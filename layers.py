"""Four analytical layers and synthesis for land-degradation detection."""
from pathlib import Path
import json
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.features import geometry_mask
from rasterio.transform import from_bounds, array_bounds
from rasterio.io import MemoryFile
from rasterio.warp import transform_bounds
import requests
from scipy.ndimage import uniform_filter, uniform_filter1d
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import csv
from datetime import datetime
import warnings

def _safe_mean(a, axis=0):
    """NaN-aware mean without warnings for all-missing slices."""
    a = np.asarray(a, dtype=np.float32)
    count = np.isfinite(a).sum(axis=axis)
    total = np.nansum(a, axis=axis)
    return np.divide(total, count, out=np.full(total.shape, np.nan, dtype=np.float32), where=count > 0)

def _safe_std(a, axis=0):
    """NaN-aware standard deviation; all-missing slices remain NaN."""
    mean = _safe_mean(a, axis=axis)
    expanded = np.expand_dims(mean, axis=axis)
    count = np.isfinite(a).sum(axis=axis)
    sq = np.nansum((a - expanded) ** 2, axis=axis)
    var = np.divide(sq, count, out=np.full(sq.shape, np.nan, dtype=np.float32), where=count > 0)
    return np.sqrt(var)

def _read(path):
    with rasterio.open(path) as src:
        a=src.read().astype("float32"); profile=src.profile.copy()
    a[a <= -9998]=np.nan
    return a,profile

def _mk_sen(y,x,chunk=20000):
    n,p=y.shape
    ii,jj=np.triu_indices(n,1)
    nv=np.isfinite(y).sum(0)
    z=np.zeros(p,dtype="float32"); slope=np.full(p,np.nan,dtype="float32")
    for s in range(0,p,chunk):
        e=min(s+chunk,p)
        d=y[jj,s:e]-y[ii,s:e]
        S=np.nansum(np.sign(d),axis=0)
        m=nv[s:e].astype(float)
        den=np.sqrt(np.maximum(m*(m-1)*(2*m+5)/18,1))
        z[s:e]=np.where(S>0,S-1,np.where(S<0,S+1,0))/den
        with np.errstate(invalid="ignore",divide="ignore"):
            slopes=d/(x[jj]-x[ii])[:,None]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                slope[s:e]=np.nanmedian(slopes,axis=0)
    return {"p":np.nan_to_num(2*(1-norm.cdf(np.abs(z))),nan=1.0),"slope":slope,"n":nv}

def run_layers(config):
    data=Path(config["data_dir"]); out=Path(config["out_dir"]); out.mkdir(parents=True,exist_ok=True)
    start,end=config["start_year"],config["end_year"]
    ndvi_path=data/"ndvi_preprocessed.tif"
    if not ndvi_path.exists(): ndvi_path=data/"ndvi.tif"
    ndvi,profile=_read(ndvi_path); rain,_=_read(data/"rain.tif"); vh,_=_read(data/"vh.tif")
    T,H,W=ndvi.shape; Y=T//12
    if T%12 or rain.shape[0]!=T: raise ValueError("NDVI and rainfall must contain matching complete monthly years.")
    years=np.arange(start,start+Y,dtype=float)
    sar=np.full(ndvi.shape,np.nan,dtype="float32")
    offset=T-vh.shape[0]
    if offset<0: raise ValueError("SAR stack cannot be longer than NDVI stack.")
    sar[offset:]=vh
    px_km2=abs(profile["transform"].a*profile["transform"].e)/1e6
    boundary=data/"boundary"/"boundary.shp"
    shp=gpd.read_file(boundary).to_crs(profile["crs"])
    inside=geometry_mask(shp.geometry,out_shape=(H,W),transform=profile["transform"],invert=True).reshape(-1)
    nd,rn,sr=(a.reshape(T,-1) for a in (ndvi,rain,sar))
    mean_ndvi=_safe_mean(nd,axis=0)
    valid=inside & np.isfinite(nd).any(0) & np.isfinite(rn).any(0) & (mean_ndvi>=config["min_mean_ndvi"])
    if valid.sum()==0: raise RuntimeError("No valid vegetated pixels; check the boundary and min_mean_ndvi.")
    nd,rn,sr=nd[:,valid],rn[:,valid],sr[:,valid]
    def grid(a,fill=0):
        z=np.full(H*W,fill,dtype="float32"); z[valid]=a
        return z.reshape(H,W)
    tr=profile["transform"]
    extent=(tr.c,tr.c+tr.a*W,tr.f+tr.e*H,tr.f)
    bg=np.where(inside.reshape(H,W),mean_ndvi.reshape(H,W),np.nan)
    finite=bg[np.isfinite(bg)]
    vlo,vhi=(np.percentile(finite,[2,98]) if finite.size else (0,1))
    def overlay(ax,cat,colors,title):
        ax.imshow(bg,cmap="Greys",extent=extent,vmin=vlo,vmax=vhi)
        ax.imshow(np.ma.masked_equal(cat,0),cmap=ListedColormap(colors),vmin=.5,vmax=len(colors)+.5,
                  extent=extent,interpolation="nearest")
        shp.boundary.plot(ax=ax,color="black",linewidth=.6)
        ax.set(xlim=extent[:2],ylim=extent[2:],xticks=[],yticks=[],title=title)
    # Annual series: only interpolate short/internal NDVI gaps in preprocessing.
    ndvi_monthly = nd.reshape(Y,12,-1)
    ndvi_month_count = np.isfinite(ndvi_monthly).sum(axis=1)
    ndvi_a = _safe_mean(ndvi_monthly,axis=1)
    ndvi_a[ndvi_month_count < 6] = np.nan
    rain_monthly = rn.reshape(Y,12,-1)
    rain_month_count = np.isfinite(rain_monthly).sum(axis=1)
    rain_a = np.nansum(rain_monthly,axis=1)
    # Annual rainfall is only comparable when all 12 monthly totals are present.
    rain_a[rain_month_count < 12] = np.nan
    sar_months=config["sar_months"]
    if sar_months is None:
        clim=_safe_mean(_safe_mean(rn.reshape(Y,12,-1),axis=2),axis=0)
        m0=int(np.argmin([sum(clim[(m+k)%12] for k in range(3)) for m in range(12)]))
        sar_months=[(m0+k)%12+1 for k in range(3)]
    season=sr.reshape(Y,12,-1)[:,np.array(sar_months)-1,:]
    # Average only where observations exist; avoid RuntimeWarning for all-NaN pixels.
    season_count=np.isfinite(season).sum(axis=1)
    season_sum=np.nansum(season,axis=1)
    sar_a=np.divide(season_sum, season_count,
                    out=np.full(season_sum.shape, np.nan, dtype="float32"),
                    where=season_count > 0)
    # Require at least two valid observations in the selected seasonal window.
    sar_a[season_count < 2] = np.nan
    r1=_mk_sen(ndvi_a,years)
    flag1=(r1["p"]<config["alpha"])&(r1["n"]>=5)&(r1["slope"]<config["ndvi_min_slope"])
    fig,ax=plt.subplots(figsize=(6,6)); overlay(ax,grid(flag1),["#d62728"],"L1 NDVI trend"); fig.savefig(out/"layer1_ndvi_trend.png",dpi=130,bbox_inches="tight"); plt.close(fig)
    # Layer 2: drought recovery lag from monthly deseasonalised NDVI.
    monthly=nd.copy()
    seasonal=_safe_mean(monthly.reshape(Y,12,-1),axis=0)
    deseason=monthly-np.tile(seasonal,(Y,1))
    sm=uniform_filter1d(deseason,size=3,axis=0,mode="nearest")
    sd=_safe_std(sm,axis=0)
    rz=(rain_a-_safe_mean(rain_a,axis=0))/(_safe_std(rain_a,axis=0)+1e-12)
    drought=rz<config["drought_z"]
    clean=sm.copy()
    for y in range(Y):
        start_m=y*12; stop_m=min(T,start_m+24)
        clean[start_m:stop_m,drought[y]]=np.nan
    lags=np.full((Y,sm.shape[1]),np.nan,dtype="float32")
    for y in range(1,Y):
        idx=np.where(drought[y]&~drought[y-1])[0]
        if not idx.size: continue
        s0=y*12
        hist=clean[max(0,s0-60):s0,idx]
        hist_count=np.isfinite(hist).sum(axis=0)
        hist_sum=np.nansum(hist,axis=0)
        base=np.divide(hist_sum, hist_count,
                       out=np.full(hist_sum.shape, np.nan, dtype="float32"),
                       where=hist_count > 0)
        window_end=min(T,s0+18)
        if window_end<=s0: continue
        local=sm[s0:window_end,idx]
        # nanargmin raises when a column has no valid values. Filter those pixels
        # before calling it, and also require enough valid historical observations.
        usable=np.isfinite(local).any(axis=0) & np.isfinite(base) & (hist_count>=6)
        if not np.any(usable):
            continue
        candidate_idx=idx[usable]
        candidate_hist=hist[:,usable]
        candidate_base=base[usable]
        candidate_local=local[:,usable]
        mi=s0+np.nanargmin(candidate_local,axis=0)
        depth=candidate_base-sm[mi,candidate_idx]
        lagvals=np.full(candidate_idx.size,37,dtype=float)
        for j,pix in enumerate(candidate_idx):
            recovery_end=min(T,mi[j]+37)
            recovered=np.where(sm[mi[j]+1:recovery_end,pix]>=candidate_base[j]-0.2*depth[j])[0]
            if recovered.size: lagvals[j]=recovered[0]+1
        ok=(mi+36<T)&(hist_count[usable]>=6)&(depth>0.5*sd[candidate_idx])
        lags[y,candidate_idx[ok]]=lagvals[ok]
    n_events=np.isfinite(lags).sum(0)
    full=np.full((Y,H*W),np.nan,dtype="float32"); full[:,valid]=lags
    full=full.reshape(Y,H,W); has=np.isfinite(full)
    smooth_lags=np.where(has,uniform_filter(np.where(has,full,0),(1,5,5),mode="constant")/
                          np.maximum(uniform_filter(has.astype("float32"),(1,5,5),mode="constant"),1e-8),np.nan)
    lags=smooth_lags.reshape(Y,-1)[:,valid]
    with np.errstate(invalid="ignore"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            lags-=np.nanmedian(lags,axis=1,keepdims=True)
    r2=_mk_sen(lags,years)
    flag2=(r2["p"]<config["lag_alpha"])&(n_events>=config["min_events"])&(r2["slope"]>config["lag_min_slope"])
    fig,ax=plt.subplots(figsize=(6,6)); overlay(ax,grid(flag2),["#9467bd"],"L2 Recovery lag"); fig.savefig(out/"layer2_recovery_lag.png",dpi=130,bbox_inches="tight"); plt.close(fig)
    # Layer 3: rainfall-adjusted NDVI residual trend (RESTREND-style linear adjustment).
    # Pixel-wise rainfall regression uses paired valid annual observations only.
    pair_ok = np.isfinite(rain_a) & np.isfinite(ndvi_a)
    pair_count = pair_ok.sum(axis=0)
    xmean = np.divide(np.where(pair_ok, rain_a, 0).sum(axis=0), pair_count,
                      out=np.full(pair_count.shape, np.nan, dtype="float32"), where=pair_count > 0)
    ymean = np.divide(np.where(pair_ok, ndvi_a, 0).sum(axis=0), pair_count,
                      out=np.full(pair_count.shape, np.nan, dtype="float32"), where=pair_count > 0)
    xm = rain_a - xmean[None, :]; ym = ndvi_a - ymean[None, :]
    cov = np.where(pair_ok, xm * ym, 0).sum(axis=0)
    var = np.where(pair_ok, xm * xm, 0).sum(axis=0)
    beta = np.divide(cov, var, out=np.full(var.shape, np.nan, dtype="float32"), where=(var > 1e-12) & (pair_count >= 5))
    resid = ym - beta[None, :] * xm
    resid[~pair_ok] = np.nan
    r3=_mk_sen(resid,years)
    flag3=(r3["p"]<config["alpha"])&(r3["n"]>=5)&(r3["slope"]<config["resid_min_slope"])
    fig,ax=plt.subplots(figsize=(6,6)); overlay(ax,grid(flag3),["#ff7f0e"],"L3 RESTREND residual"); fig.savefig(out/"layer3_restrend.png",dpi=130,bbox_inches="tight"); plt.close(fig)
    # Layer 4: annual dry-season Sentinel-1 VH structure trend.
    have=np.isfinite(sar_a).any(1)
    if have.sum()<5: raise RuntimeError("Fewer than five years of usable Sentinel-1 observations.")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        common=np.nanmedian(sar_a[have],axis=1)
    r4=_mk_sen(sar_a[have]-common[:,None],years[have])
    flag4=(r4["p"]<config["alpha"])&(r4["n"]>=5)&(r4["slope"]<config["sar_min_slope"])
    r4n=_mk_sen(ndvi_a[have],years[have])
    sar_masked=flag4&~((r4n["p"]<config["alpha"])&(r4n["slope"]<0))
    fig,ax=plt.subplots(figsize=(6,6)); overlay(ax,grid(flag4),["#1f77b4"],"L4 SAR structure"); fig.savefig(out/"layer4_sar_structure.png",dpi=130,bbox_inches="tight"); plt.close(fig)
    # Explicitly name the four binary evidence layers for downstream comparison.
    flags = {"L1_NDVI_trend": flag1, "L2_recovery_lag": flag2,
             "L3_RESTREND": flag3, "L4_SAR_structure": flag4}
    confidence=sum(f.astype("uint8") for f in flags.values())
    degraded=confidence>0
    # Alternative-method comparison on the same valid pixels. These are
    # detection-footprint comparisons, not accuracy scores (no reference labels).
    baseline_ndvi = flag1
    climate_adjusted = flag1 | flag3
    optical_resilience = flag1 | flag2 | flag3
    four_layer = flag1 | flag2 | flag3 | flag4
    comparison_methods = {
        "NDVI_only": baseline_ndvi,
        "NDVI_plus_RESTREND": climate_adjusted,
        "Optical_plus_resilience": optical_resilience,
        "Four_layer_fusion": four_layer,
    }
    reference_metrics = None
    reference_path = config.get("reference_raster")
    if reference_path:
        ref_path = Path(reference_path)
        if not ref_path.exists():
            raise FileNotFoundError(f"Configured reference_raster does not exist: {ref_path}")
        with rasterio.open(ref_path) as ref_src:
            if (ref_src.height, ref_src.width) != (H, W) or ref_src.transform != profile["transform"] or ref_src.crs != profile["crs"]:
                raise ValueError("Reference raster must be aligned to the analysis grid (same CRS, transform, width and height). Reproject/resample it before validation.")
            ref = ref_src.read(1).astype("float32")
            ref_nodata = ref_src.nodata
        ref_flat = ref.reshape(-1)[valid]
        ref_ok = np.isfinite(ref_flat)
        if ref_nodata is not None:
            ref_ok &= ref_flat != ref_nodata
        if not np.isin(ref_flat[ref_ok], [0, 1]).all():
            raise ValueError("Reference raster must encode 1=degraded and 0=not degraded; NoData is ignored.")
        if ref_ok.sum() == 0 or np.unique(ref_flat[ref_ok]).size < 2:
            raise ValueError("Reference raster needs valid examples of both degraded (1) and not-degraded (0) pixels.")
        reference_metrics = {}
        for method_name, prediction in comparison_methods.items():
            pred = prediction[ref_ok].astype(bool); truth = ref_flat[ref_ok].astype(bool)
            tp = int(np.sum(pred & truth)); fp = int(np.sum(pred & ~truth))
            fn = int(np.sum(~pred & truth)); tn = int(np.sum(~pred & ~truth))
            precision = tp / max(tp + fp, 1); recall = tp / max(tp + fn, 1)
            specificity = tn / max(tn + fp, 1); f1 = 2 * precision * recall / max(precision + recall, 1e-12)
            reference_metrics[method_name] = {
                "n_reference_pixels": int(ref_ok.sum()), "TP": tp, "FP": fp, "FN": fn, "TN": tn,
                "accuracy": (tp + tn) / max(tp + fp + fn + tn, 1),
                "precision": precision, "recall": recall, "specificity": specificity, "F1": f1,
                "IoU": tp / max(tp + fp + fn, 1), "balanced_accuracy": (recall + specificity) / 2,
            }
    driver=np.zeros(confidence.shape,dtype="uint8")
    # Mutually exclusive explanatory labels with explicit precedence:
    # recovery-lag evidence first, then residual/SAR evidence, then NDVI-only.
    driver[flag1 & ~(flag2 | flag3 | flag4)] = 1
    driver[(flag3 | flag4) & ~flag2] = 2
    driver[flag2] = 3
    missed=degraded&~flag1
    km2=lambda m: round(float(np.sum(m))*px_km2,2)
    summary={
      "study_period":[int(start),int(end)],"valid_pixels":int(valid.sum()),"valid_area_km2":km2(np.ones(valid.sum(),bool)),
      "degraded_area_km2":km2(degraded),"degraded_pct_of_valid":round(100*float(degraded.mean()),2),
      "detected_by_baseline_km2":km2(flag1),"missed_by_baseline_km2":km2(missed),
      "missed_pct_of_degraded":round(100*missed.sum()/max(degraded.sum(),1),2),
      "by_layer_km2":{"L1 NDVI trend":km2(flag1),"L2 Recovery lag":km2(flag2),"L3 RESTREND residual":km2(flag3),"L4 SAR structure":km2(flag4)},
      "by_driver_km2":{"NDVI-only signal":km2(driver==1),"Rainfall-adjusted or SAR signal (heuristic)":km2(driver==2),"Recovery-lag signal (priority heuristic)":km2(driver==3)},
      "high_confidence_km2":km2(confidence>=2),"sar_down_ndvi_stable_km2":km2(sar_masked),
      "diagnostics":{"median_ndvi_slope_per_year":float(np.nanmedian(r1["slope"])) if np.isfinite(r1["slope"]).any() else None,
                     "median_recovery_lag_slope_months_per_year":float(np.nanmedian(r2["slope"])) if np.isfinite(r2["slope"]).any() else None,
                     "median_restrend_slope_per_year":float(np.nanmedian(r3["slope"])) if np.isfinite(r3["slope"]).any() else None,
                     "median_sar_relative_slope_per_year":float(np.nanmedian(r4["slope"])) if np.isfinite(r4["slope"]).any() else None,
                     "sar_years_used":int(have.sum()), "sar_months":sar_months,
                     "pct_pixels_with_enough_drought_events":round(100*float((n_events>=config["min_events"]).mean()),1),
                     "pixels_with_valid_ndvi":int(np.isfinite(ndvi_a).any(0).sum()),
                     "pixels_with_valid_rainfall":int(np.isfinite(rain_a).any(0).sum()),
                     "pixels_with_valid_sar":int(np.isfinite(sar_a).any(0).sum())},
      "method_comparison": {
          "NDVI_only": {"area_km2": km2(baseline_ndvi), "pixels_flagged": int(baseline_ndvi.sum()),
                         "interpretation": "Simple greenness-trend baseline; can miss resilience, rainfall-adjusted and radar-only signals."},
          "NDVI_plus_RESTREND": {"area_km2": km2(climate_adjusted), "pixels_flagged": int(climate_adjusted.sum()),
                         "interpretation": "Adds a simple rainfall-adjusted residual trend to the NDVI baseline."},
          "Optical_plus_resilience": {"area_km2": km2(optical_resilience), "pixels_flagged": int(optical_resilience.sum()),
                         "interpretation": "Adds recovery-lag evidence but remains dependent on optical greenness."},
          "Four_layer_fusion": {"area_km2": km2(four_layer), "pixels_flagged": int(four_layer.sum()),
                         "interpretation": "Combines greenness, recovery, rainfall-adjusted residuals and SAR structural-change evidence."},
          "additional_area_vs_NDVI_km2": km2(four_layer & ~baseline_ndvi),
          "additional_pixels_vs_NDVI": int((four_layer & ~baseline_ndvi).sum()),
          "reference_validation_available": reference_metrics is not None,
          "reference_metrics": reference_metrics,
          "comparison_warning": "Larger mapped area is not proof of higher accuracy. Use independent reference labels or field validation to compare precision, recall, F1 and IoU."
      }
    }
    # Export per-layer continuous trend surfaces and per-layer diagnostics so
    # analysts can interpret effect size, significance and sample support.
    slope_arrays = [r1["slope"], r2["slope"], r3["slope"], r4["slope"]]
    p_arrays = [r1["p"], r2["p"], r3["p"], r4["p"]]
    slope_names = ["L1_NDVI_SenSlope_per_year", "L2_recovery_lag_slope_per_year",
                   "L3_RESTREND_residual_slope_per_year", "L4_SAR_relative_slope_per_year"]
    rows_out = []
    for name, flag, rr, unit in zip(slope_names, flags.values(),
                                     [r1, r2, r3, r4],
                                     ["NDVI/year", "months/year", "NDVI residual/year", "relative VH units/year"]):
        finite_slope = rr["slope"][np.isfinite(rr["slope"])]
        rows_out.append({"layer": name, "flagged_pixels": int(flag.sum()), "flagged_area_km2": km2(flag),
                         "valid_slope_pixels": int(np.isfinite(rr["slope"]).sum()),
                         "median_slope": float(np.median(finite_slope)) if finite_slope.size else None,
                         "slope_unit": unit, "p_threshold": float(config.get("lag_alpha", config["alpha"]) if "recovery" in name else config["alpha"])})
    with (out / "layer_diagnostics.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_out[0].keys())); writer.writeheader(); writer.writerows(rows_out)
    comp_rows = []
    for key in ("NDVI_only", "NDVI_plus_RESTREND", "Optical_plus_resilience", "Four_layer_fusion"):
        val = summary["method_comparison"][key]
        metrics = (reference_metrics or {}).get(key, {})
        comp_rows.append({"method": key, "area_km2": val["area_km2"], "pixels_flagged": val["pixels_flagged"],
                              "precision": metrics.get("precision"), "recall": metrics.get("recall"),
                              "F1": metrics.get("F1"), "IoU": metrics.get("IoU"),
                              "balanced_accuracy": metrics.get("balanced_accuracy"),
                              "reference_pixels": metrics.get("n_reference_pixels"), "interpretation": val["interpretation"]})
    with (out / "method_comparison.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["method", "area_km2", "pixels_flagged", "precision", "recall", "F1", "IoU", "balanced_accuracy", "reference_pixels", "interpretation"]); writer.writeheader(); writer.writerows(comp_rows)
    # Four-band trend product: useful in GIS for comparing effect size, not just flags.
    slope_profile = {**profile, "count": 4, "dtype": "float32", "nodata": -9999, "compress": "lzw"}
    with rasterio.open(out / "layer_trend_slopes.tif", "w", **slope_profile) as dst:
        for i, arr in enumerate(slope_arrays, start=1):
            dst.write(np.nan_to_num(grid(arr, np.nan), nan=-9999).astype("float32"), i)
            dst.set_band_description(i, slope_names[i-1])
    # Human-readable report: explicitly distinguish observations from hypotheses.
    report = ["# Land Degradation Detection — Interpretive Results", "",
              f"**Study period:** {start}–{end}", f"**Valid analysis pixels:** {int(valid.sum()):,}",
              f"**Estimated valid area:** {summary['valid_area_km2']:.2f} km²",
              f"**Pixels flagged by at least one layer:** {int(degraded.sum()):,} ({summary['degraded_pct_of_valid']:.2f}% of valid pixels)", "",
              "## How to interpret the four layers", "",
              "- **L1 NDVI trend:** persistent decline in annual vegetation greenness. This is the baseline and is relatively easy to explain, but it can miss structural changes and resilience loss.",
              "- **L2 recovery lag:** increasing time for vegetation anomalies to recover after rainfall-defined dry years. Treat as a resilience indicator; event count and missing observations matter.",
              "- **L3 RESTREND-style residual:** declining NDVI after a simple linear rainfall adjustment. It suggests a decline not explained by this model's rainfall relationship, not proof of a human cause.",
              "- **L4 SAR structure:** declining relative Sentinel-1 VH backscatter in the selected seasonal window. This can add structural evidence but can also respond to soil moisture and surface conditions.", "",
              "## Comparison with alternative approaches", "",
              "| Method | Flagged area (km²) | What it adds |", "|---|---:|---|"]
    for key in ("NDVI_only", "NDVI_plus_RESTREND", "Optical_plus_resilience", "Four_layer_fusion"):
        val = summary["method_comparison"][key]
        report.append(f"| {key.replace('_',' ')} | {val['area_km2']:.2f} | {val['interpretation']} |")
    report += ["", f"Four-layer area not flagged by NDVI-only baseline: **{summary['method_comparison']['additional_area_vs_NDVI_km2']:.2f} km²**.",
               "This is a difference in mapped detection footprint, not an accuracy gain. A larger footprint can include true positives and false positives.", ""]
    if reference_metrics:
        report += ["Independent reference labels were supplied, so classification metrics were calculated on aligned valid pixels:", "",
                   "| Method | Precision | Recall | F1 | IoU | Balanced accuracy | Reference pixels |", "|---|---:|---:|---:|---:|---:|---:|"]
        for key, metrics in reference_metrics.items():
            report.append(f"| {key.replace('_',' ')} | {metrics['precision']:.3f} | {metrics['recall']:.3f} | {metrics['F1']:.3f} | {metrics['IoU']:.3f} | {metrics['balanced_accuracy']:.3f} | {metrics['n_reference_pixels']} |")
        report += ["", "Use these scores to compare alternatives on the same reference pixels. Check class balance and spatial/temporal independence before making a superiority claim.", ""]
    else:
        report += ["No independent reference raster was configured. Therefore this run compares mapped area and detection footprint only; it cannot establish which method is more accurate. To enable accuracy comparison, set `reference_raster` in `run_pipeline.py` to an aligned raster coded 1=degraded, 0=not degraded, and NoData elsewhere.", ""]
    report += ["## Confidence score", "",
               "The confidence raster counts how many of the four layer flags are positive (0–4). It is an evidence-agreement score, not a probability. Layer outputs can be correlated because several depend on NDVI or rainfall.", "",
               "## Limitations and recommended validation", "",
               "1. Thresholds are configurable and should be sensitivity-tested; results may change with alpha, slope cut-offs, minimum event counts and the study period.",
               "2. The recovery-lag method uses simplified drought and recovery definitions. Inspect representative time series and compare events with known drought records.",
               "3. The RESTREND implementation is a simple rainfall-adjusted linear residual trend; temperature, land cover, irrigation and other drivers are not modelled. Annual rainfall is used only when all 12 monthly totals are present.",
               "4. Sentinel-1 VH change is not uniquely attributable to canopy degradation; soil moisture, roughness and acquisition geometry can affect it.",
               "5. If independent reference labels are unavailable, do not claim higher accuracy than alternatives. Use field data, expert-labelled samples or an independent degradation product for validation.",
               "6. Compare alternatives on the same valid pixels and period. Report precision, recall, F1, IoU and area-adjusted accuracy only after reference labels exist.", "",
               "## Output guide", "",
               "- `flagged_region.png`: side-by-side layer flags, synthesis, heuristic driver, NDVI-only omissions and area comparison.",
               "- `confidence.tif`: 0–4 layer-agreement count.", "- `driver.tif`: heuristic category only; not causal attribution.",
               "- `layer_trend_slopes.tif`: continuous slope/effect-size surfaces for all four layers.",
               "- `layer_diagnostics.csv`: pixel counts, mapped area, valid slope support and median slope per layer.",
               "- `method_comparison.csv`: comparable detection footprints and, when reference labels are configured, precision/recall/F1/IoU metrics.",
               "- `summary.json`: machine-readable metrics and diagnostics.", "",
               f"Generated: {datetime.now().isoformat(timespec='seconds')}" ]
    (out / "interpretation_report.md").write_text("\n".join(report), encoding="utf-8")
    fig,ax=plt.subplots(figsize=(8,4.5))
    method_labels=[r["method"].replace("_","\n") for r in comp_rows]
    method_areas=[r["area_km2"] for r in comp_rows]
    bars=ax.bar(method_labels, method_areas)
    ax.set_ylabel("Flagged area (km²)"); ax.set_title("Alternative-method detection footprint comparison")
    ax.tick_params(axis="x", labelsize=8); ax.grid(axis="y", alpha=.25)
    for bar, value in zip(bars, method_areas): ax.text(bar.get_x()+bar.get_width()/2, bar.get_height(), f"{value:.1f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout(); fig.savefig(out / "method_comparison.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    # Data coverage chart makes it clear where each input supports the analysis.
    coverage_labels = ["NDVI", "Rainfall", "SAR"]
    coverage_values = [float(np.isfinite(ndvi_a).any(0).mean()*100),
                       float(np.isfinite(rain_a).any(0).mean()*100),
                       float(np.isfinite(sar_a).any(0).mean()*100)]
    fig, axcov = plt.subplots(figsize=(6.5, 4))
    bars = axcov.bar(coverage_labels, coverage_values)
    axcov.set_ylim(0, 100); axcov.set_ylabel("Valid analysis pixels (%)")
    axcov.set_title("Input data coverage within eligible pixels")
    axcov.grid(axis="y", alpha=.25)
    for bar, value in zip(bars, coverage_values):
        axcov.text(bar.get_x()+bar.get_width()/2, bar.get_height(), f"{value:.1f}%", ha="center", va="bottom")
    fig.tight_layout(); fig.savefig(out / "data_coverage.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    # Use the exact Leaflet basemap screenshot captured in the browser at the
    # user's chosen zoom and viewport. Reproject the analysis flags onto the
    # screenshot pixel grid (not the other way around) to preserve image detail.
    leaflet_rgb = None
    leaflet_extent = None
    leaflet_transform = None
    leaflet_inside = None
    leaflet_shape = None
    capture_path = config.get("leaflet_capture_path")
    capture_bounds = config.get("leaflet_bounds_3857")
    if capture_path and Path(capture_path).exists() and capture_bounds and len(capture_bounds) == 4:
        try:
            from PIL import Image
            with Image.open(capture_path) as captured:
                image = np.asarray(captured.convert("RGB"), dtype=np.uint8)
            ch, cw = image.shape[:2]
            west, south, east, north = map(float, capture_bounds)
            if east <= west or north <= south:
                raise ValueError("Invalid Leaflet screenshot bounds.")
            leaflet_transform = from_bounds(west, south, east, north, cw, ch)
            leaflet_extent = (west, east, south, north)
            leaflet_shape = (ch, cw)
            # Transform the original user polygon to the Leaflet map's Web Mercator CRS.
            shp_3857 = shp.to_crs("EPSG:3857")
            leaflet_inside_full = geometry_mask(shp_3857.geometry, out_shape=leaflet_shape,
                                                 transform=leaflet_transform, invert=True)
            yy, xx = np.where(leaflet_inside_full)
            if not len(xx):
                raise ValueError("The drawn pentagon does not overlap the captured Leaflet map extent.")
            # Crop to the pixel bounding box of the drawn polygon, then retain a precise polygon alpha mask.
            x0, x1 = int(xx.min()), int(xx.max()) + 1
            y0, y1 = int(yy.min()), int(yy.max()) + 1
            from affine import Affine
            leaflet_transform = leaflet_transform * Affine.translation(x0, y0)
            leaflet_rgb = image[y0:y1, x0:x1].astype("float32") / 255.0
            leaflet_shape = (y1-y0, x1-x0)
            leaflet_inside = leaflet_inside_full[y0:y1, x0:x1]
            leaflet_extent = tuple(array_bounds(leaflet_shape[0], leaflet_shape[1], leaflet_transform))
            # array_bounds returns west,south,east,north; matplotlib expects left,right,bottom,top.
            west, south, east, north = leaflet_extent
            leaflet_extent = (west, east, south, north)
            print(f"Using cropped Leaflet pentagon screenshot ({leaflet_shape[1]}x{leaflet_shape[0]}) at the user's drawn map zoom.")
        except Exception as capture_error:
            print(f"Warning: Leaflet screenshot could not be used: {capture_error}")

    def original_flag_overlay(ax, flag, color, title):
        if leaflet_rgb is not None:
            # Transparent outside the exact pentagon. Preserve screenshot pixels at native capture resolution.
            base_rgba = np.dstack([leaflet_rgb, leaflet_inside.astype("float32")])
            ax.imshow(base_rgba, extent=leaflet_extent, interpolation="nearest", origin="upper")
            # Reproject the categorical detection flags onto the screenshot grid with nearest-neighbour.
            flag_src = grid(flag).astype("uint8")
            flag_dst = np.zeros(leaflet_shape, dtype="uint8")
            reproject(source=flag_src, destination=flag_dst,
                      src_transform=profile["transform"], src_crs=profile["crs"],
                      dst_transform=leaflet_transform, dst_crs="EPSG:3857",
                      resampling=Resampling.nearest, src_nodata=0, dst_nodata=0)
            flag_dst[~leaflet_inside] = 0
            rgba = np.zeros((*leaflet_shape, 4), dtype="float32")
            from matplotlib.colors import to_rgba
            rgba[flag_dst > 0] = to_rgba(color, alpha=0.62)
            ax.imshow(rgba, extent=leaflet_extent, interpolation="nearest", origin="upper")
            shp_3857.boundary.plot(ax=ax, color="#d9f477", linewidth=1.5)
            ax.set(xlim=leaflet_extent[:2], ylim=leaflet_extent[2:], xticks=[], yticks=[], title=title)
        else:
            # Fallback for API-created runs that did not originate from the Leaflet UI.
            ax.imshow(bg, cmap="Greys", extent=extent, interpolation="nearest")
            mask = (grid(flag) > 0) & inside.reshape(H, W)
            rgba = np.zeros((H, W, 4), dtype="float32")
            from matplotlib.colors import to_rgba
            rgba[mask] = to_rgba(color, alpha=0.58)
            ax.imshow(rgba, extent=extent, interpolation="nearest")
            shp.boundary.plot(ax=ax, color="#d9f477", linewidth=1.2)
            ax.set(xlim=extent[:2], ylim=extent[2:], xticks=[], yticks=[], title=title)

    fig,ax=plt.subplots(2,4,figsize=(18,9))
    layer_specs = [(flag1,"#d62728","L1 NDVI trend"),(flag2,"#9467bd","L2 Recovery lag"),
                   (flag3,"#ff7f0e","L3 RESTREND"),(flag4,"#1f77b4","L4 SAR structure")]
    for a,(f,c,nm) in zip(ax[0], layer_specs):
        original_flag_overlay(a, f, c, nm + " · Leaflet satellite basemap")
        # Save each layer as a standalone original-image detection output too.
        one_fig, one_ax = plt.subplots(figsize=(8, 8))
        original_flag_overlay(one_ax, f, c, nm + " · Leaflet satellite basemap")
        one_fig.tight_layout()
        one_fig.savefig(out / (nm.lower().replace(" ", "_").replace("·", "").strip("_") + "_original_flagged.png"),
                        dpi=160, bbox_inches="tight")
        plt.close(one_fig)
    overlay(ax[1,0],grid(confidence),["#ffe066","#ffa94d","#f03e3e","#7a0010"],"Degradation confidence")
    overlay(ax[1,1],grid(driver),["#e6a117","#c0392b","#6c5ce7"],"Likely driver (heuristic)")
    overlay(ax[1,2],grid(missed),["#8e44ad"],"Missed by NDVI-only baseline")
    base,miss=summary["detected_by_baseline_km2"],summary["missed_by_baseline_km2"]
    layered_area=summary["degraded_area_km2"]
    ax[1,3].bar(["NDVI only","Four-layer"],[base,layered_area])
    ax[1,3].set(ylabel="Area (km²)",title=f"Additional area: {miss:.2f} km²")
    fig.tight_layout(); fig.savefig(out/"flagged_region.png",dpi=140,bbox_inches="tight"); plt.close(fig)
    out_profile={**profile,"count":1,"dtype":"float32","nodata":-9999,"compress":"lzw"}
    for name,a in {"confidence":confidence,"driver":driver,"missed":missed,"layer1_ndvi":flag1,"layer2_recovery":flag2,"layer3_restrend":flag3,"layer4_sar":flag4}.items():
        raster_values = np.full(H*W, np.nan, dtype="float32")
        raster_values[valid] = np.asarray(a, dtype="float32")
        with rasterio.open(out/f"{name}.tif","w",**out_profile) as dst:
            dst.write(np.nan_to_num(raster_values.reshape(H,W),nan=-9999).astype("float32"),1)
    (out/"summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
