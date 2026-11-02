"""Automatic regional reference retrieval from the public Trends.Earth SDG 15.3.1 COG.

The source is a global cloud-optimized GeoTIFF hosted by Zenodo. GDAL/rasterio
uses HTTP range requests to read only the selected study-area window.
"""
from __future__ import annotations

import json
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import rasterio
import geopandas as gpd
from rasterio.windows import from_bounds
from rasterio.warp import reproject, Resampling
from rasterio.enums import Resampling as RioResampling

ZENODO_RECORD = "https://zenodo.org/records/17514520"
COG_URL = (
    "https://zenodo.org/records/17514520/files/"
    "TrendsEarth_SDG15.3.1_2000-2023_Trends.Earth.tif"
)
SOURCE_BAND = 14  # 2023 status relative to 2000-2015 baseline
SOURCE_NAME = "Trends.Earth SDG Indicator 15.3.1 (default LPD)"


def create_regional_reference(boundary_file: str, target_grid: str, output_path: str) -> dict:
    """Read Band 14 for the AOI and align it to target_grid.

    Trends.Earth status codes: 1-3 = degradation; 4-7 = stability or
    improvement; -32768 = NoData. This is an external benchmark/proxy, not
    field-verified ground truth.
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    boundary_gdf = gpd.read_file(boundary_file)
    if boundary_gdf.empty or boundary_gdf.crs is None:
        raise ValueError("Study boundary is empty or has no CRS; cannot query reference data.")
    bounds = boundary_gdf.total_bounds
    boundary_crs = boundary_gdf.crs

    with rasterio.open(target_grid) as grid:
        target_profile = grid.profile.copy()
        dst_crs = grid.crs
        dst_transform = grid.transform
        dst_height, dst_width = grid.height, grid.width
        dst_bounds = grid.bounds

    # /vsicurl/ enables GDAL HTTP range reads; the source is a COG, so only
    # the requested regional window and required overview levels are fetched.
    vsi_url = "/vsicurl/" + COG_URL
    try:
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif", CPL_VSIL_CURL_USE_HEAD="NO"):
            with rasterio.open(vsi_url) as src:
                if src.count < SOURCE_BAND:
                    raise RuntimeError(f"Reference COG has {src.count} bands; expected band {SOURCE_BAND}.")
                if src.crs is None:
                    raise RuntimeError("Reference COG has no CRS metadata.")
                from rasterio.warp import transform_bounds
                ref_bounds = transform_bounds(boundary_crs, src.crs, *bounds, densify_pts=21)
                win = from_bounds(*ref_bounds, transform=src.transform)
                win = win.round_offsets().round_lengths()
                # Clamp to the global raster extent.
                col0 = max(0, int(win.col_off)); row0 = max(0, int(win.row_off))
                col1 = min(src.width, int(win.col_off + win.width)); row1 = min(src.height, int(win.row_off + win.height))
                if col1 <= col0 or row1 <= row0:
                    raise ValueError("The selected study region is outside the Trends.Earth reference coverage.")
                window = rasterio.windows.Window(col0, row0, col1-col0, row1-row0)
                raw = src.read(SOURCE_BAND, window=window, masked=True)
                src_transform = src.window_transform(window)
                src_crs = src.crs
                source_nodata = src.nodata
    except Exception as exc:
        raise RuntimeError(
            "Could not read the regional Trends.Earth COG. Check internet access, "
            "Zenodo availability and GDAL's HTTP range support. Original error: " + str(exc)
        ) from exc

    raw_arr = np.asarray(raw.filled(-32768), dtype=np.int16)
    # Keep only documented SDG status classes. Unknown codes remain NoData.
    regional = np.full(raw_arr.shape, 255, dtype=np.uint8)
    regional[np.isin(raw_arr, [1, 2, 3])] = 1
    regional[np.isin(raw_arr, [4, 5, 6, 7])] = 0
    if not np.any(regional != 255):
        raise ValueError("Trends.Earth has no valid reference pixels in this study region.")

    aligned = np.full((dst_height, dst_width), 255, dtype=np.uint8)
    reproject(
        source=regional,
        destination=aligned,
        src_transform=src_transform,
        src_crs=src_crs,
        src_nodata=255,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        dst_nodata=255,
        resampling=Resampling.nearest,
    )
    target_profile.update(count=1, dtype="uint8", nodata=255, compress="deflate")
    with rasterio.open(out, "w", **target_profile) as dst:
        dst.write(aligned, 1)

    valid = aligned != 255
    degraded = int(np.count_nonzero(aligned == 1))
    stable = int(np.count_nonzero(aligned == 0))
    metadata = {
        "source": SOURCE_NAME,
        "source_url": COG_URL,
        "record_url": ZENODO_RECORD,
        "license": "CC BY 4.0",
        "source_band": SOURCE_BAND,
        "source_period": "2023 status compared with 2000-2015 baseline",
        "source_codes": {"1": "persistent degradation", "2": "recent degradation", "3": "baseline degradation", "4": "stability", "5": "baseline improvement", "6": "recent improvement", "7": "persistent improvement", "-32768": "NoData"},
        "binary_mapping": {"1": "degraded (source status 1-3)", "0": "not degraded (source status 4-7)", "255": "NoData"},
        "reference_type": "external benchmark proxy; not field-verified ground truth",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "valid_pixels": int(valid.sum()),
        "degraded_pixels": degraded,
        "not_degraded_pixels": stable,
        "nodata_pixels": int((~valid).sum()),
        "valid_coverage_percent": round(100 * float(valid.mean()), 2),
        "raster": str(out),
    }
    out.with_suffix(".json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
