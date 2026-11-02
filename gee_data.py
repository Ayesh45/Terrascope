"""Google Earth Engine retrieval for the four-layer degradation pipeline."""
from pathlib import Path
import time
import numpy as np
import geopandas as gpd
import rasterio
import json
from datetime import datetime
from rasterio.crs import CRS

def retrieve_data(config):
    import ee, requests
    boundary_file = config["boundary_file"]
    data_dir = Path(config["data_dir"]); data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "boundary").mkdir(parents=True, exist_ok=True)
    g = gpd.read_file(boundary_file)
    if g.crs is None:
        g = g.set_crs(4326)
    g = g[g.geom_type.isin(["Polygon", "MultiPolygon"])].to_crs(4326).dissolve()
    if g.empty:
        raise ValueError("Boundary file contains no polygon geometry.")
    g.to_file(data_dir / "boundary" / "boundary.shp", driver="ESRI Shapefile")
    geom = g.geometry.iloc[0].simplify(0.001)
    crs = g.estimate_utm_crs().to_string()
    x0, y0, x1, y1 = gpd.GeoSeries([geom], crs=4326).to_crs(crs).total_bounds
    scale = int(config["scale"])
    rows, cols = int(np.ceil((y1-y0)/scale))+1, int(np.ceil((x1-x0)/scale))+1
    start, end, sar_start = config["start_year"], config["end_year"], config["sar_start_year"]
    ee_project = config["ee_project"]
    try:
        ee.Initialize(project=ee_project)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=ee_project)
    region = ee.Geometry(geom.__geo_interface__)

    def safe(col, band, reducer):
        empty = ee.Image.constant(0).rename(band).updateMask(ee.Image.constant(0))
        return ee.Image(ee.Algorithms.If(col.size().gt(0), col.select(band).reduce(reducer).rename(band), empty))
    def ndvi_month(s,e):
        c = ee.ImageCollection("MODIS/061/MOD13Q1").filterDate(s,e).filterBounds(region)
        c = c.map(lambda im: im.select("NDVI").updateMask(im.select("SummaryQA").lte(1)))
        return safe(c, "NDVI", ee.Reducer.mean()).multiply(0.0001)
    def rain_month(s,e):
        return safe(ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY").filterDate(s,e).filterBounds(region),
                    "precipitation", ee.Reducer.sum())
    def vh_month(s,e):
        c = (ee.ImageCollection("COPERNICUS/S1_GRD").filterDate(s,e).filterBounds(region)
             .filter(ee.Filter.eq("instrumentMode","IW"))
             .filter(ee.Filter.listContains("transmitterReceiverPolarisation","VH"))
             .filter(ee.Filter.eq("orbitProperties_pass","DESCENDING")))
        return safe(c, "VH", ee.Reducer.median())
    from preprocessing import mask_s2_clouds
    def s2_month(s,e):
        c = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterDate(s,e).filterBounds(region)
             .map(mask_s2_clouds))
        # NDVI after SCL + QA60 cloud/cirrus masking. Empty months remain masked.
        empty = ee.Image.constant(0).rename("NDVI").updateMask(ee.Image.constant(0))
        def calc(im):
            return im.normalizedDifference(["nir","red"]).rename("NDVI").copyProperties(im,["system:time_start"])
        c = c.map(calc)
        return ee.Image(ee.Algorithms.If(c.size().gt(0), c.median().rename("NDVI"), empty))
    sources = {
        "ndvi": (ndvi_month, start, 36),
        "rain": (rain_month, start, 36),
        "vh": (vh_month, sar_start, 12),
        "s2_ndvi": (s2_month, max(2017, start), 12),
    }
    if rows*cols*4 > 40_000_000:
        raise RuntimeError("One band exceeds Earth Engine download limit; increase SCALE.")
    for name, (fn, first_year, max_bands) in sources.items():
        target = data_dir / f"{name}.tif"
        expected_count = 12 * (end - first_year + 1)
        if target.exists():
            with rasterio.open(target) as src:
                l, b, r, t = src.bounds
                bounds_ok = max(abs(l-x0), abs(b-y0), abs(r-x1), abs(t-y1)) <= 3*scale
                if (src.count == expected_count and src.crs == CRS.from_string(crs)
                        and bounds_ok and src.width == cols and src.height == rows):
                    print(f"Using existing {target}")
                    continue
                raise RuntimeError(f"{target} does not match current area/grid/period. Move it or remove it and rerun.")
        months = [(y,m) for y in range(first_year,end+1) for m in range(1,13)]
        n = max(1, min(max_bands, 40_000_000 // (rows*cols*4)))
        chunks=[]
        for k in range(0,len(months),n):
            chunk_path = data_dir / f".{name}_chunk_{k//n:03d}.tif"
            if not chunk_path.exists():
                print(f"{name}: months {k+1}-{min(k+n,len(months))}/{len(months)}")
                images=[]
                for y,m in months[k:k+n]:
                    s=ee.Date.fromYMD(y,m,1); e=s.advance(1,"month")
                    images.append(fn(s,e).rename(f"{y}_{m:02d}"))
                img=ee.Image.cat(images)
                if name in ("vh","s2_ndvi"):
                    img=(img.setDefaultProjection(ee.Projection(crs).atScale(10 if name=="vh" else 10))
                         .reduceResolution(ee.Reducer.mean(), maxPixels=1024))
                url=img.toFloat().unmask(-9999).getDownloadURL(
                    {"region":region,"crs":crs,"scale":scale,"format":"GEO_TIFF"})
                for attempt in range(3):
                    resp=requests.get(url,timeout=900)
                    if resp.content[:4] in (b"II*\x00",b"MM\x00*"): break
                    print("Retrying GEE download:",resp.text[:160]); time.sleep(4*(attempt+1))
                else:
                    raise RuntimeError(f"Earth Engine download failed for {name}: {resp.text[:300]}")
                chunk_path.write_bytes(resp.content)
            chunks.append(chunk_path)
        arrays=[]
        for f in chunks:
            with rasterio.open(f) as src:
                arrays.append(src.read().astype("float32")); profile=src.profile.copy()
        stack=np.concatenate(arrays,axis=0)
        if stack.shape[0] != expected_count:
            raise RuntimeError(f"{name}: expected {expected_count} monthly bands, got {stack.shape[0]}")
        profile.update(count=stack.shape[0],dtype="float32",nodata=-9999,compress="lzw")
        with rasterio.open(target,"w",**profile) as dst: dst.write(stack)
        for f in chunks: f.unlink(missing_ok=True)
        print("Saved",target,stack.shape)
    # Download an unenhanced Sentinel-2 true-colour composite for displaying
    # detections over the original satellite appearance (not the CLAHE preview).
    rgb_target = data_dir / "s2_rgb_original.tif"
    if not rgb_target.exists():
        rgb_start = ee.Date.fromYMD(max(2017, start), 1, 1)
        rgb_end = ee.Date.fromYMD(end + 1, 1, 1)
        rgb_collection = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterDate(rgb_start, rgb_end).filterBounds(region)
            .map(mask_s2_clouds))
        rgb_count = rgb_collection.size().getInfo()
        if rgb_count == 0:
            raise RuntimeError("No cloud-masked Sentinel-2 scenes are available for an original RGB image in this region and period.")
        rgb = (rgb_collection.median().select(["red", "green", "blue"], ["R", "G", "B"])
               .clip(region).toFloat())
        rgb_url = rgb.unmask(-9999).getDownloadURL({
            "region": region, "crs": crs, "scale": scale, "format": "GEO_TIFF"
        })
        rgb_resp = requests.get(rgb_url, timeout=900)
        if rgb_resp.status_code != 200 or rgb_resp.content[:4] not in (b"II*\x00", b"MM\x00*"):
            raise RuntimeError(f"Sentinel-2 original RGB download failed: {rgb_resp.text[:300]}")
        rgb_target.write_bytes(rgb_resp.content)
        print("Saved original Sentinel-2 true-colour image:", rgb_target)

    # Assert all analysis stacks share a grid; S2 is optional recent imagery but same grid.
    grids=set()
    for name in ("ndvi","rain","vh"):
        with rasterio.open(data_dir/f"{name}.tif") as src:
            grids.add((src.transform,src.width,src.height,src.crs))
    if len(grids)!=1:
        raise RuntimeError("Downloaded analysis stacks do not share the same raster grid.")
    # Keep a compact data dictionary beside the actual downloaded stacks. This
    # makes source, temporal coverage, units and analytical purpose inspectable.
    catalog = {
      "created_at": datetime.now().isoformat(timespec="seconds"),
      "study_area": str(boundary_file), "crs": crs,
      "grid": {"rows": rows, "columns": cols, "requested_scale_m": scale,
               "bounds_projected": [float(x) for x in (x0,y0,x1,y1)]},
      "datasets": {
        "ndvi.tif": {"collection": "MODIS/061/MOD13Q1", "bands": "monthly NDVI composites from product NDVI band, SummaryQA <= 1",
          "period": [start, end], "units": "unitless NDVI, scaled by 0.0001", "used_for": ["Layer 1 annual greenness trend", "Layer 2 recovery anomalies", "Layer 3 rainfall-adjusted residual trend"],
          "limitations": "MODIS is a coarser-resolution long-term baseline; quality filtering does not eliminate every atmospheric or compositing artefact."},
        "rain.tif": {"collection": "UCSB-CHG/CHIRPS/DAILY", "bands": "monthly sum of daily precipitation", "period": [start, end],
          "units": "mm/month", "used_for": ["drought-year detection in Layer 2", "rainfall adjustment in Layer 3", "dry-season window selection for Layer 4 when sar_months is null"],
          "limitations": "Rainfall is a gridded estimate and may not represent local microclimates or irrigation."},
        "vh.tif": {"collection": "COPERNICUS/S1_GRD", "bands": "monthly median VH backscatter, descending IW acquisitions", "period": [sar_start, end],
          "units": "dB backscatter; seasonal trend is relative VH dB/year after annual spatial-median adjustment", "used_for": ["Layer 4 relative seasonal structural-change trend"],
          "limitations": "Backscatter also responds to moisture, surface roughness, geometry and acquisition conditions."},
        "s2_ndvi.tif": {"collection": "COPERNICUS/S2_SR_HARMONIZED", "bands": "monthly median NDVI after SCL and QA60 cloud/cirrus masking", "period": [max(2017, start), end],
          "units": "unitless NDVI", "used_for": ["recent high-resolution visual preview and image-quality context"],
          "limitations": "This supplementary stack is not substituted for the MODIS long-term series in the four analytical layers."}
      },
      "missing_value": -9999,
      "important_note": "NoData represents unavailable observations, not zero vegetation, zero rainfall or no degradation."
    }
    (data_dir / "data_catalog.json").write_text(json.dumps(catalog, indent=2), encoding="utf-8")
    print("Saved data dictionary:", data_dir / "data_catalog.json")
    return {"crs":crs,"bounds":[float(x) for x in (x0,y0,x1,y1)],"rows":rows,"cols":cols,
            "catalog": str(data_dir / "data_catalog.json")}
