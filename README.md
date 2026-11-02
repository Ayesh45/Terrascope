# TerraScope — Land Degradation Explorer

A local web application for drawing a five-vertex study polygon, downloading its ESRI Shapefile bundle, and running the existing four-layer land-degradation pipeline.

## Features

- Satellite basemap with polygon drawing/editing.
- Requires exactly five vertices before submission.
- Exports the study boundary as GeoJSON and a ZIP of Shapefile sidecars (`.shp`, `.shx`, `.dbf`, `.prj`, etc.).
- Runs the included Earth Engine retrieval, preprocessing, four-layer analysis and alternative-method comparison.
- Displays processing progress/logs and provides artifact download links.
- Each job has its own folder under `runs/` to avoid mixing separate study areas.

## Setup (Windows PowerShell)

Use Python 3.10–3.12 in a virtual environment. GeoPandas/Rasterio are native geospatial packages; if pip has trouble installing them on Windows, use a Conda environment.

```powershell
cd land_degradation_webapp
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
$env:EE_PROJECT="your-earth-engine-enabled-cloud-project-id"
python -c "import ee; ee.Authenticate()"
python app.py
```

Then open `http://127.0.0.1:8000`.

### Google Earth Engine prerequisites

1. Have a Google account with Earth Engine access and an enabled Google Cloud project.
2. Set `EE_PROJECT` to that project ID before starting the server.
3. Run `python -c "import ee; ee.Authenticate()"` once in the same environment. Follow the browser authentication flow.
4. Ensure the Earth Engine project is registered/configured for use with Earth Engine.

The application does not include API keys or credentials. Do not place credentials in the frontend. On a deployed server, configure service-account credentials or another approved non-interactive authentication method; never commit a private key.

## Run the analysis

1. Open the website and zoom to your area.
2. Select the polygon drawing tool from the map toolbar.
3. Click five points around the study area. Finish the polygon by clicking the first point (Leaflet Draw may also finish on double-click; ensure it has five distinct vertices).
4. Select years and pixel scale.
5. Click **Run degradation analysis**.
6. Download the shapefile ZIP immediately, then wait for processing. When complete, use the artifact links to download the reports and rasters.

## Important interpretation notes

- The shapefile is generated before satellite processing starts.
- Earth Engine data retrieval may take several minutes or longer depending on area and period. Large polygons can hit Earth Engine request/download limits; use a smaller area or larger pixel size if necessary.
- The processing worker is single-job-at-a-time for this prototype. Job status is held in memory; restarting the server clears the dashboard's job registry, though files remain in `runs/`.
- The confidence raster represents the number of analytical layer flags, not a calibrated probability. Driver categories are heuristic indicators, not definitive causal attribution.
- Alternative-method mapped-area comparisons do not prove that one method is more accurate. Independent reference labels are required for accuracy metrics.
- This research prototype has no user accounts, multi-user authorization, or production job queue. Do not expose it publicly without adding authentication, rate limits, storage controls, and deployment-specific Earth Engine credentials.

## Main files

- `app.py`: FastAPI API and background job orchestration.
- `templates/index.html`: web page structure.
- `static/app.js`: map, drawing, job polling and result links.
- `static/style.css`: responsive interface.
- `gee_data.py`, `preprocessing.py`, `layers.py`: research processing pipeline.


## Independent reference raster and validation

After a landscape analysis completes, the website displays an **Independent validation** panel. Create labelled polygons in QGIS (or another GIS) and upload them as a GeoJSON file or a ZIP containing a shapefile and its `.shp`, `.shx`, `.dbf`, and `.prj` files. The polygon attribute named in **Class attribute name** must contain either `1` / `degraded` for degraded land or `0` / `not degraded` / `non-degraded` for non-degraded land. The default attribute name is `class`.

The website reprojects the polygons to the analysis CRS and rasterizes them onto the exact grid of `confidence.tif`. Areas outside labelled polygons are stored as NoData (255). It then reruns the local layer-analysis stage using the already downloaded/preprocessed data; it does **not** repeat Earth Engine retrieval. The generated `reference_upload/reference.tif` is available in the downloads, and `method_comparison.csv`, `summary.json`, and the interpretation report are updated with reference-based metrics for all four methods.

**Reference quality matters:** labels must be independently interpreted from imagery, field observations, or a suitable external inventory. Do not create the reference labels from the model's own output. Include both classes, use polygons only, and ensure the labels represent the same time period as the analysis. The rasterization uses pixel-centre inclusion (`all_touched=False`) and the model's exact grid. Inspect the raster before reporting metrics.

Supported reference uploads: `.geojson`, `.json`, or `.zip` shapefile bundle; maximum upload size 50 MB.

## Automatic regional external reference benchmark

After the four-layer screening run, TerraScope attempts to retrieve the study-area window from the public Trends.Earth SDG Indicator 15.3.1 global Cloud-Optimized GeoTIFF (Zenodo record `17514520`). The implementation uses GDAL/rasterio `/vsicurl/` HTTP range reads so it should read only the regional COG window rather than downloading the complete multi-gigabyte global file.

- Source: `https://zenodo.org/records/17514520`
- COG: `TrendsEarth_SDG15.3.1_2000-2023_Trends.Earth.tif`
- License: CC BY 4.0
- Band 14: 2023 status compared with the 2000–2015 baseline.
- Source status codes 1–3 are mapped to `1` (degradation); codes 4–7 are mapped to `0` (not degraded); `-32768` and unknown codes remain NoData (`255`).
- The regional raster is reprojected/aligned to the exact analysis grid using nearest-neighbour resampling.
- If retrieval fails, the main four-layer analysis still completes and reports the reference source as unavailable. If the overlap contains only one class, the raster is retained but binary classification metrics are not calculated.

**Scientific caveat:** Trends.Earth is an independently published SDG 15.3.1 benchmark, not field-verified ground truth. Its land-productivity component has methodological overlap with TerraScope's NDVI-based layers. Treat resulting accuracy statistics as external benchmark agreement, disclose the 2023 reference period, and do not describe them as field-validated accuracy. ESA WorldCover and SoilGrids are not used as binary degradation labels: they provide land-cover and soil context, respectively.

The regional reference file and its provenance metadata are saved under each run's `reference_auto/` directory. The web interface displays its coverage and evaluation status and provides a download link.

### Leaflet screenshot-based flagged-region visuals

When an analysis is started from the web interface, TerraScope captures the currently visible Esri World Imagery tiles from the Leaflet map at the user's selected zoom and viewport. The browser sends this PNG and its Web Mercator bounds with the study polygon. The analysis then crops the capture to the polygon's pixel bounds, masks pixels outside the exact polygon, and reprojects each analytical flag layer onto the screenshot pixel grid using nearest-neighbour resampling. The four individual flagged-region PNGs and the combined `flagged_region.png` therefore use the actual Leaflet map capture rather than a separately fetched basemap or Sentinel-2/preprocessed image. This is a visual-only change; analytical rasters and calculations are unchanged.

The Leaflet tile layer uses cross-origin image loading for browser screenshot capture. If a browser extension, network policy, or tile-server CORS response prevents capture, the UI reports the capture error rather than silently substituting a different satellite image.
