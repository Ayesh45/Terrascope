from __future__ import annotations
import io, json, os, re, shutil, threading, traceback, uuid, zipfile, tempfile, base64
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime, timezone
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import geopandas as gpd
from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from shapely.geometry import shape, mapping, Polygon, MultiPolygon
from rasterio.features import rasterize
import rasterio

from gee_data import retrieve_data
from preprocessing import preprocess_downloaded_data
from layers import run_layers
from reference_data import create_regional_reference

ROOT = Path(__file__).resolve().parent
WORK = ROOT / "runs"
WORK.mkdir(exist_ok=True)
EE_PROJECT = os.getenv("EE_PROJECT", "")
executor = ThreadPoolExecutor(max_workers=1)
lock = threading.Lock()
jobs: dict[str, dict] = {}

app = FastAPI(title="Land Degradation Explorer", version="1.0.0")
app.mount("/static", StaticFiles(directory=ROOT / "static"), name="static")

class RunRequest(BaseModel):
    geometry: dict
    start_year: int = Field(default=2001, ge=2001, le=2025)
    end_year: int = Field(default=2025, ge=2002, le=2025)
    scale: int = Field(default=250, ge=100, le=1000)
    ee_project: str | None = None
    leaflet_capture: str | None = None
    leaflet_bounds_3857: list[float] | None = None

@app.get("/", response_class=HTMLResponse)
def home():
    return FileResponse(ROOT / "templates" / "index.html")

@app.get("/api/health")
def health():
    return {"ok": True, "earth_engine_project_configured": bool(EE_PROJECT)}

@app.post("/api/jobs")
def create_job(req: RunRequest):
    folder = None  # initialized before validation/capture handling
    if req.end_year <= req.start_year:
        raise HTTPException(400, "End year must be later than start year.")
    try:
        geom = shape(req.geometry)
    except Exception as e:
        raise HTTPException(400, f"Invalid polygon geometry: {e}")
    if geom.geom_type != "Polygon" or geom.is_empty or not geom.is_valid:
        raise HTTPException(400, "Draw a valid, non-self-intersecting polygon.")
    # Polygon exterior contains the closing coordinate as its final repeated point.
    vertices = len(geom.exterior.coords) - 1
    if vertices != 5:
        raise HTTPException(400, f"Please draw exactly five vertices (a pentagon). Current vertices: {vertices}.")
    minx, miny, maxx, maxy = geom.bounds
    if minx < -180 or maxx > 180 or miny < -85 or maxy > 85:
        raise HTTPException(400, "Polygon coordinates must be valid longitude/latitude values.")
    ee_project = (req.ee_project or EE_PROJECT).strip()
    if not ee_project:
        raise HTTPException(400, "Google Earth Engine project is not configured. Set EE_PROJECT in the environment before starting the website.")
    # Create the job workspace before saving the browser capture.
    job_id = uuid.uuid4().hex[:12]
    folder = WORK / job_id
    data_dir, out_dir = folder / "data", folder / "outputs"
    data_dir.mkdir(parents=True); out_dir.mkdir(parents=True)

    # Preserve the exact Leaflet basemap capture made in the browser at the drawn map zoom.
    leaflet_capture_path = None
    if req.leaflet_capture and req.leaflet_bounds_3857 and len(req.leaflet_bounds_3857) == 4:
        try:
            encoded = req.leaflet_capture.split(",", 1)[-1]
            image_bytes = base64.b64decode(encoded, validate=True)
            if len(image_bytes) > 18 * 1024 * 1024:
                raise ValueError("Leaflet map capture exceeds 18 MB.")
            from PIL import Image
            with Image.open(io.BytesIO(image_bytes)) as captured:
                if captured.format != "PNG" or captured.width < 100 or captured.height < 100:
                    raise ValueError("Leaflet capture must be a valid PNG image.")
                captured.verify()
            leaflet_capture_path = folder / "leaflet_map_capture.png"
            leaflet_capture_path.write_bytes(image_bytes)
        except Exception as capture_error:
            shutil.rmtree(folder, ignore_errors=True)
            raise HTTPException(400, f"Could not read the Leaflet map capture: {capture_error}")
    boundary_geojson = folder / "study_area.geojson"
    boundary_geojson.write_text(json.dumps({"type":"FeatureCollection","features":[{"type":"Feature","properties":{"name":"study_area"},"geometry":mapping(geom)}]}, indent=2))
    # Generate an ESRI shapefile bundle immediately so the user can download it before analysis finishes.
    gdf = gpd.GeoDataFrame({"name":["study_area"]}, geometry=[geom], crs="EPSG:4326")
    shp_dir = folder / "shapefile"; shp_dir.mkdir()
    shp_path = shp_dir / "study_area.shp"
    gdf.to_file(shp_path, driver="ESRI Shapefile", encoding="UTF-8")
    shp_zip = folder / "study_area_shapefile.zip"
    with zipfile.ZipFile(shp_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in shp_dir.iterdir(): zf.write(p, p.name)
    config = {
        "boundary_file": str(boundary_geojson),
        "ee_project": ee_project,
        "start_year": req.start_year, "end_year": req.end_year,
        "sar_start_year": max(2017, req.start_year), "scale": req.scale,
        "data_dir": str(data_dir), "out_dir": str(out_dir),
        "min_mean_ndvi": 0.10, "alpha": 0.05, "ndvi_min_slope": -0.0005,
        "drought_z": -1.0, "min_events": 4, "lag_alpha": 0.10,
        "lag_min_slope": 0.1, "resid_min_slope": -0.001, "sar_min_slope": -0.05,
        "sar_months": None, "reference_raster": None,
        "leaflet_capture_path": str(leaflet_capture_path) if leaflet_capture_path else None,
        "leaflet_bounds_3857": req.leaflet_bounds_3857,
    }
    with lock:
        jobs[job_id] = {"id":job_id,"status":"queued","stage":"Queued","progress":0,"message":"Waiting to start","created_at":datetime.now(timezone.utc).isoformat(),"error":None,"folder":str(folder),"config":config}
    executor.submit(run_job, job_id)
    return {"job_id":job_id,"status":"queued","shapefile_url":f"/api/jobs/{job_id}/download/study_area_shapefile.zip","status_url":f"/api/jobs/{job_id}"}

def run_job(job_id: str):
    with lock:
        job = jobs[job_id]; job["status"]="running"; job["stage"]="Starting"; job["progress"]=3
    folder = Path(job["folder"]); config = job["config"]
    log_path = folder / "process.log"
    class LogWriter(io.TextIOBase):
        def write(self, s):
            if not s: return 0
            with open(log_path, "a", encoding="utf-8") as f: f.write(s)
            low=s.lower()
            with lock:
                j=jobs[job_id]
                if "step 1" in low or "months " in low or "retrieving" in low:
                    j["stage"]="Retrieving satellite and rainfall data"; j["progress"]=max(j["progress"],10)
                elif "step 2" in low or "preprocess" in low:
                    j["stage"]="Preprocessing data"; j["progress"]=max(j["progress"],65)
                elif "step 3" in low or "layer" in low or "analysis" in low:
                    j["stage"]="Calculating degradation layers and comparisons"; j["progress"]=max(j["progress"],78)
                elif "completed" in low:
                    j["stage"]="Complete"; j["progress"]=100
            return len(s)
    stream=LogWriter()
    try:
        with redirect_stdout(stream), redirect_stderr(stream):
            print("Step 1/3: Retrieving monthly satellite and rainfall data")
            retrieval=retrieve_data(config)
            with lock: jobs[job_id]["progress"]=65; jobs[job_id]["stage"]="Preprocessing data"
            print("Step 2/3: Preprocessing downloaded data")
            qa=preprocess_downloaded_data(config["data_dir"], config["out_dir"])
            with lock: jobs[job_id]["progress"]=78; jobs[job_id]["stage"]="Calculating layers and comparing methods"
            print("Step 3/4: Running four-layer analysis and method comparison")
            summary=run_layers(config)
            automatic_reference = None
            try:
                with lock:
                    jobs[job_id]["stage"] = "Retrieving regional external reference data"
                    jobs[job_id]["progress"] = 88
                print("Step 4/4: Retrieving regional Trends.Earth SDG 15.3.1 reference data")
                ref_path = folder / "reference_auto" / "trends_earth_reference.tif"
                automatic_reference = create_regional_reference(
                    config["boundary_file"], str(Path(config["out_dir"]) / "confidence.tif"), str(ref_path)
                )
                config["reference_raster"] = str(ref_path)
                try:
                    summary = run_layers(config)
                    automatic_reference["evaluation_status"] = "metrics_calculated"
                    automatic_reference["evaluation_note"] = "Metrics compare TerraScope outputs with the external Trends.Earth SDG 15.3.1 benchmark. They are not field-validated accuracy scores."
                except ValueError as validation_error:
                    # A tiny AOI or one-class reference footprint cannot support
                    # binary classification metrics. Keep the reference raster
                    # and complete the run without claiming accuracy metrics.
                    config["reference_raster"] = None
                    summary = run_layers(config)
                    automatic_reference["evaluation_status"] = "reference_downloaded_metrics_unavailable"
                    automatic_reference["evaluation_note"] = str(validation_error)
                except Exception:
                    config["reference_raster"] = None
                    raise
            except Exception as reference_error:
                # Reference data is supplementary: an outage or unavailable
                # region must not invalidate the core four-layer analysis.
                config["reference_raster"] = None
                automatic_reference = {
                    "source": "Trends.Earth SDG Indicator 15.3.1",
                    "evaluation_status": "unavailable",
                    "error": str(reference_error),
                    "note": "Core TerraScope analysis completed without external reference evaluation."
                }
                print("Automatic reference retrieval unavailable:", reference_error)
            (folder/"run_metadata.json").write_text(json.dumps({"retrieval":retrieval,"preprocessing_qa":qa,"summary":summary,"automatic_reference":automatic_reference,"config":{k:v for k,v in config.items() if k not in ("ee_project",)}}, indent=2, default=str))
        with lock:
            jobs[job_id].update(status="completed",stage="Complete",progress=100,message="Analysis completed.")
    except Exception as e:
        with open(log_path,"a",encoding="utf-8") as f: f.write("\nERROR\n"+traceback.format_exc())
        with lock:
            jobs[job_id].update(status="failed",stage="Failed",message=str(e),error=traceback.format_exc())

@app.post("/api/jobs/{job_id}/reference")
async def upload_reference(job_id: str, file: UploadFile = File(...), class_field: str = Form("class")):
    """Convert independently labelled polygons into a reference raster on the model grid, then evaluate all methods."""
    with lock:
        job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found in this server session.")
    if job.get("status") != "completed":
        raise HTTPException(409, "Wait until the analysis has completed before adding reference labels.")
    filename = (file.filename or "").lower()
    if not (filename.endswith(".geojson") or filename.endswith(".json") or filename.endswith(".zip")):
        raise HTTPException(400, "Upload a GeoJSON (.geojson/.json) or zipped shapefile (.zip).")
    folder = Path(job["folder"]); out_dir = Path(job["config"]["out_dir"])
    grid_candidates = [out_dir / "confidence.tif"]
    grid_path = next((p for p in grid_candidates if p.exists()), None)
    if grid_path is None:
        raise HTTPException(500, "The confidence raster is missing; cannot determine the analysis grid.")
    ref_dir = folder / "reference_upload"; ref_dir.mkdir(exist_ok=True)
    payload = await file.read()
    if len(payload) > 50 * 1024 * 1024:
        raise HTTPException(413, "Reference upload exceeds the 50 MB limit.")
    source_path = ref_dir / (Path(file.filename or "reference.geojson").name)
    source_path.write_bytes(payload)
    try:
        if filename.endswith(".zip"):
            extract_dir = ref_dir / "shapefile"; extract_dir.mkdir(exist_ok=True)
            with zipfile.ZipFile(source_path) as zf:
                members = [m for m in zf.namelist() if not m.startswith("__MACOSX/") and ".." not in Path(m).parts]
                zf.extractall(extract_dir, members)
            shp_files = list(extract_dir.rglob("*.shp"))
            if not shp_files:
                raise ValueError("ZIP must contain a .shp file and its .shx/.dbf/.prj companions.")
            vector_path = shp_files[0]
        else:
            vector_path = source_path
        gdf = gpd.read_file(vector_path)
        if gdf.empty or gdf.crs is None:
            raise ValueError("Reference polygons must be non-empty and have a defined CRS.")
        if class_field not in gdf.columns:
            raise ValueError(f"Class field '{class_field}' was not found. Available fields: {', '.join(map(str, gdf.columns))}")
        gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
        if gdf.empty or not gdf.geom_type.isin(["Polygon", "MultiPolygon"]).all():
            raise ValueError("The reference file must contain polygon or multipolygon features only.")
        values = gdf[class_field].astype(str).str.strip().str.lower().map({"0":0, "not degraded":0, "non-degraded":0, "non_degraded":0, "1":1, "degraded":1})
        if values.isna().any():
            raise ValueError("Class values must be 0/1 or 'not degraded'/'degraded'.")
        gdf["_class"] = values.astype("uint8")
        with rasterio.open(grid_path) as grid:
            profile = grid.profile.copy(); transform = grid.transform; crs = grid.crs
            height, width = grid.height, grid.width
        projected = gdf.to_crs(crs)
        shapes = [(geom, int(value)) for geom, value in zip(projected.geometry, projected["_class"]) if geom is not None and not geom.is_empty]
        ref_array = rasterize(shapes, out_shape=(height, width), transform=transform, fill=255, dtype="uint8", all_touched=False)
        present = set(int(v) for v in set(ref_array.ravel()) if v in (0,1))
        if present != {0,1}:
            raise ValueError("The rasterized reference must cover both classes (0=not degraded and 1=degraded) on valid analysis pixels. Check polygon coverage and labels.")
        ref_path = ref_dir / "reference.tif"
        profile.update(count=1, dtype="uint8", nodata=255, compress="deflate")
        with rasterio.open(ref_path, "w", **profile) as dst:
            dst.write(ref_array, 1)
        config = dict(job["config"]); config["reference_raster"] = str(ref_path)
        # Re-run only the local analysis stage; Earth Engine downloads and preprocessing are reused.
        summary = run_layers(config)
        metadata_path = folder / "run_metadata.json"
        metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
        metadata["summary"] = summary
        metadata["config"]["reference_raster"] = "reference_upload/reference.tif"
        metadata_path.write_text(json.dumps(metadata, indent=2, default=str))
        with lock:
            jobs[job_id]["config"] = config
            jobs[job_id]["reference_ready"] = True
        return {"ok": True, "message": "Reference raster created and all four methods evaluated.", "reference_raster": f"/api/jobs/{job_id}/download/reference_upload/reference.tif", "summary": summary.get("method_comparison", {}).get("reference_metrics")}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(400, f"Could not create/evaluate reference raster: {e}")


@app.get("/api/jobs/{job_id}/reference-info")
def reference_info(job_id: str):
    with lock:
        job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    metadata_path = Path(job["folder"]) / "run_metadata.json"
    if not metadata_path.exists():
        return {"status": "pending"}
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    info = metadata.get("automatic_reference") or {"evaluation_status": "not_available"}
    ref_file = Path(job["folder"]) / "reference_auto" / "trends_earth_reference.tif"
    if ref_file.exists():
        info["download_url"] = f"/api/jobs/{job_id}/download/reference_auto/trends_earth_reference.tif"
    return info

@app.get("/api/jobs/{job_id}")
def get_job(job_id: str):
    with lock:
        if job_id not in jobs: raise HTTPException(404,"Job not found in this server session.")
        j=jobs[job_id].copy()
    folder=Path(j["folder"])
    artifacts=[]
    # Traverse the job folder exactly once. Scanning folder, outputs, and data
    # separately repeats every nested output because folder.rglob already descends into them.
    seen_paths=set()
    if folder.exists():
        for p in folder.rglob("*"):
            if not p.is_file() or p.name == "process.log" or p.name.startswith("."):
                continue
            rel=p.relative_to(folder).as_posix()
            if rel in seen_paths:
                continue
            seen_paths.add(rel)
            artifacts.append({"name":p.name,"path":rel,"url":f"/api/jobs/{job_id}/download/{rel}"})
    artifacts.sort(key=lambda item: (item["name"].lower(), item["path"].lower()))
    j.pop("folder",None); j.pop("config",None); j.pop("error",None)
    j["artifacts"]=artifacts
    if (folder/"process.log").exists():
        j["log"]=(folder/"process.log").read_text(encoding="utf-8",errors="replace")[-12000:]
    return j

@app.get("/api/jobs/{job_id}/download/{file_path:path}")
def download(job_id: str, file_path: str):
    with lock:
        job=jobs.get(job_id)
    if not job: raise HTTPException(404,"Job not found")
    base=Path(job["folder"]).resolve(); target=(base/file_path).resolve()
    if base not in target.parents and target != base: raise HTTPException(400,"Invalid path")
    if not target.exists() or not target.is_file(): raise HTTPException(404,"File not ready or not found")
    return FileResponse(target, filename=target.name)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app:app", host="127.0.0.1", port=int(os.getenv("PORT","8000")), reload=False)
