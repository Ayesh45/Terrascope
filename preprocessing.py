"""Preprocessing for satellite land-degradation analysis.

Cloud/shadow masking is performed in Earth Engine before Sentinel-2 compositing.
Local preprocessing applies conservative temporal despiking, short-gap interpolation,
and CLAHE only to display previews (not to scientific index values).
"""
from pathlib import Path
import numpy as np
import rasterio
import json
from datetime import datetime

def mask_s2_clouds(image):
    """Mask Sentinel-2 SR Harmonized clouds/cirrus/shadows using SCL and QA60."""
    scl = image.select("SCL")
    # Exclude saturated, cloud shadow, medium/high cloud, cirrus and snow/ice.
    scl_ok = scl.neq(1).And(scl.neq(3)).And(scl.neq(8)).And(
        scl.neq(9)).And(scl.neq(10)).And(scl.neq(11))
    qa = image.select("QA60")
    cloud_bit = 1 << 10
    cirrus_bit = 1 << 11
    qa_ok = qa.bitwiseAnd(cloud_bit).eq(0).And(qa.bitwiseAnd(cirrus_bit).eq(0))
    return image.updateMask(scl_ok).updateMask(qa_ok).select(
        ["B2", "B3", "B4", "B8"], ["blue", "green", "red", "nir"]
    ).divide(10000).copyProperties(image, ["system:time_start"])

def enhance_for_display(array):
    """Apply CLAHE to a single 2-D image for visualization only; requires OpenCV."""
    try:
        import cv2
    except ImportError as exc:
        raise ImportError("Install opencv-python to use CLAHE previews: pip install opencv-python") from exc
    a = np.asarray(array, dtype=np.float32)
    valid = np.isfinite(a)
    if not valid.any():
        return a
    lo, hi = np.nanpercentile(a[valid], [2, 98])
    if hi <= lo:
        return a
    scaled = np.clip((a - lo) / (hi - lo), 0, 1)
    u8 = np.nan_to_num(scaled * 255, nan=0).astype("uint8")
    result = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8)).apply(u8).astype("float32") / 255
    result[~valid] = np.nan
    return result

def preprocess_stack(stack, max_spike=0.35):
    """Remove isolated monthly NDVI spikes and interpolate only internal short gaps.

    Values remain on the original NDVI scale. Long leading/trailing gaps stay NaN.
    """
    arr = np.asarray(stack, dtype=np.float32).copy()
    arr[arr <= -9998] = np.nan
    if arr.ndim != 3:
        raise ValueError("Expected raster stack with shape (time, rows, columns)")
    # Replace isolated temporal outliers only when both neighboring observations agree.
    if arr.shape[0] >= 3:
        prev, cur, nxt = arr[:-2], arr[1:-1], arr[2:]
        isolated = np.isfinite(prev) & np.isfinite(cur) & np.isfinite(nxt)
        isolated &= (np.abs(cur - prev) > max_spike) & (np.abs(cur - nxt) > max_spike)
        cur[isolated] = (prev[isolated] + nxt[isolated]) / 2
        arr[1:-1] = cur
    t = np.arange(arr.shape[0])
    flat = arr.reshape(arr.shape[0], -1)
    for k in range(flat.shape[1]):
        ok = np.isfinite(flat[:, k])
        if ok.sum() >= 2:
            missing = ~ok
            # Fill internal gaps only, and only gaps up to 3 monthly observations.
            for idx in np.where(missing)[0]:
                left = np.where(ok & (t < idx))[0]
                right = np.where(ok & (t > idx))[0]
                if left.size and right.size:
                    l, r = left[-1], right[0]
                    if r - l <= 4:
                        flat[idx, k] = np.interp(idx, [l, r], [flat[l, k], flat[r, k]])
                        ok[idx] = True
    return flat.reshape(arr.shape)

def preprocess_downloaded_data(data_dir, out_dir):
    """Validate data stacks and create a CLAHE-enhanced preview when OpenCV is installed."""
    data_dir, out_dir = Path(data_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ndvi_path = data_dir / "ndvi.tif"
    if not ndvi_path.exists():
        return {"status": "skipped", "reason": "ndvi.tif not downloaded yet"}
    with rasterio.open(ndvi_path) as src:
        stack = src.read().astype("float32")
        profile = src.profile.copy()
    raw = stack.copy()
    raw[raw <= -9998] = np.nan
    cleaned = preprocess_stack(stack)
    raw_valid = np.isfinite(raw)
    clean_valid = np.isfinite(cleaned)
    qa = {
      "created_at": datetime.now().isoformat(timespec="seconds"),
      "input": str(ndvi_path), "output": str(data_dir / "ndvi_preprocessed.tif"),
      "shape_time_rows_cols": list(cleaned.shape),
      "input_valid_observations": int(raw_valid.sum()),
      "output_valid_observations": int(clean_valid.sum()),
      "input_missing_percent": round(100 * float((~raw_valid).mean()), 3),
      "output_missing_percent": round(100 * float((~clean_valid).mean()), 3),
      "newly_filled_observations": int((~raw_valid & clean_valid).sum()),
      "spike_threshold": 0.35,
      "interpolation_rule": "Only internal gaps bounded by valid observations and spanning at most four monthly indices are interpolated.",
      "enhancement_rule": "CLAHE is applied only to preview images; it does not alter analytical NDVI values.",
      "warning": "Spike correction and interpolation are heuristic. Validate them against representative time series before scientific interpretation."
    }
    profile.update(count=cleaned.shape[0], dtype="float32", nodata=-9999, compress="lzw")
    clean_path = data_dir / "ndvi_preprocessed.tif"
    with rasterio.open(clean_path, "w", **profile) as dst:
        dst.write(np.where(np.isfinite(cleaned), cleaned, -9999).astype("float32"))
    count = np.isfinite(cleaned).sum(axis=0)
    total = np.nansum(cleaned, axis=0)
    preview = np.divide(total, count, out=np.full(total.shape, np.nan, dtype=np.float32), where=count > 0)
    previews = []
    try:
        import matplotlib.pyplot as plt
        enhanced = enhance_for_display(preview)
        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(enhanced, cmap="YlGn", vmin=0, vmax=1)
        ax.set_title("CLAHE-enhanced mean MODIS NDVI (display only)")
        ax.axis("off")
        fig.colorbar(im, ax=ax, shrink=0.75)
        fig.tight_layout()
        modis_preview = out_dir / "preprocessing_preview.png"
        fig.savefig(modis_preview, dpi=140, bbox_inches="tight")
        plt.close(fig)
        previews.append(str(modis_preview))
        # Sentinel-2 has already been cloud/shadow masked in GEE; enhance only for display.
        s2_path = data_dir / "s2_ndvi.tif"
        if s2_path.exists():
            with rasterio.open(s2_path) as src:
                s2 = src.read().astype("float32")
            s2[s2 <= -9998] = np.nan
            s2_count = np.isfinite(s2).sum(axis=0)
            s2_total = np.nansum(s2, axis=0)
            s2_mean = np.divide(s2_total, s2_count, out=np.full(s2_total.shape, np.nan, dtype=np.float32), where=s2_count > 0)
            s2_enhanced = enhance_for_display(s2_mean)
            fig, ax = plt.subplots(figsize=(7, 6))
            im = ax.imshow(s2_enhanced, cmap="YlGn", vmin=0, vmax=1)
            ax.set_title("Cloud-masked Sentinel-2 NDVI + CLAHE (display only)")
            ax.axis("off")
            fig.colorbar(im, ax=ax, shrink=0.75)
            fig.tight_layout()
            s2_preview = out_dir / "s2_cloudmasked_clahe_preview.png"
            fig.savefig(s2_preview, dpi=140, bbox_inches="tight")
            plt.close(fig)
            previews.append(str(s2_preview))
    except ImportError:
        qa["preview_note"] = "OpenCV/matplotlib unavailable; numerical preprocessing completed but preview generation was skipped."
    qa_path = out_dir / "preprocessing_qa.json"
    qa["previews"] = previews
    qa["qa_report"] = str(qa_path)
    qa_path.write_text(json.dumps(qa, indent=2), encoding="utf-8")
    return {"status": "ok", "cleaned_ndvi": str(clean_path), "previews": previews, "qa_report": str(qa_path)}
