const map = L.map('map', {zoomControl:true}).setView([10.35,76.95], 7);
L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', {maxZoom:19, crossOrigin:'anonymous', attribution:'Tiles © Esri'}).addTo(map);
L.control.scale({imperial:false}).addTo(map);
let drawnItems = new L.FeatureGroup(); map.addLayer(drawnItems);
let currentLayer = null, currentJob = null, pollHandle = null, lastArtifactSignature = '';
const runBtn=document.getElementById('runBtn'), clearBtn=document.getElementById('clearBtn');
const startSel=document.getElementById('startYear'), endSel=document.getElementById('endYear');
for(let y=2001;y<=2025;y++){const o=document.createElement('option');o.value=y;o.textContent=y;if(y===2001)o.selected=true;startSel.appendChild(o)}
for(let y=2002;y<=2025;y++){const o=document.createElement('option');o.value=y;o.textContent=y;if(y===2025)o.selected=true;endSel.appendChild(o)}
const drawControl = new L.Control.Draw({position:'topleft',draw:{polyline:false,rectangle:false,circle:false,circlemarker:false,marker:false,polygon:{allowIntersection:false,showArea:true,shapeOptions:{color:'#d9f477',weight:3,fillColor:'#d9f477',fillOpacity:.18}}},edit:{featureGroup:drawnItems,remove:true}});map.addControl(drawControl);
map.on('mousemove',e=>document.getElementById('mapCoords').textContent=`${e.latlng.lat.toFixed(5)}°, ${e.latlng.lng.toFixed(5)}°`);
map.on(L.Draw.Event.CREATED,e=>{drawnItems.clearLayers();currentLayer=e.layer;drawnItems.addLayer(currentLayer);updateBoundary()});
map.on(L.Draw.Event.EDITED,updateBoundary);map.on(L.Draw.Event.DELETED,()=>{currentLayer=null;updateBoundary()});
function getVertices(){if(!currentLayer)return 0;const latlngs=currentLayer.getLatLngs()[0]||[];return latlngs.length>1&&latlngs[0].equals(latlngs[latlngs.length-1])?latlngs.length-1:latlngs.length}
function updateBoundary(){const n=getVertices();document.getElementById('vertexCount').textContent=`${n} / 5`;document.getElementById('boundaryStatus').textContent=!currentLayer?'Not drawn':n===5?'Pentagon ready':`${n} vertices`;runBtn.disabled=!(currentLayer&&n===5);document.getElementById('inlineError').hidden=true}
clearBtn.addEventListener('click',()=>{drawnItems.clearLayers();currentLayer=null;updateBoundary()});
runBtn.addEventListener('click',async()=>{
  const n=getVertices();if(n!==5){showError(`Draw exactly five vertices. Your polygon has ${n}.`);return}
  if(Number(endSel.value)<=Number(startSel.value)){showError('End year must be later than start year.');return}
  const feature=currentLayer.toGeoJSON();
  runBtn.disabled=true;runBtn.querySelector('span').textContent='Capturing map…';
  try{
    const capture = await captureLeafletBasemap();
    runBtn.querySelector('span').textContent='Submitting…';
    const res=await fetch('/api/jobs',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({geometry:feature.geometry,start_year:Number(startSel.value),end_year:Number(endSel.value),scale:Number(document.getElementById('scale').value),leaflet_capture:capture.dataUrl,leaflet_bounds_3857:capture.bounds3857})});
    const data=await res.json();if(!res.ok)throw new Error(data.detail||'Could not start the analysis.');
    currentJob=data.job_id;lastArtifactSignature='';document.getElementById('jobPanel').hidden=false;document.getElementById('jobPanel').scrollIntoView({behavior:'smooth',block:'start'});document.getElementById('shapefileLink').href=data.shapefile_url;document.getElementById('resultsArea').innerHTML='';document.getElementById('metricsSection').hidden=true;
    if(pollHandle)clearInterval(pollHandle);await poll();pollHandle=setInterval(poll,2500);
  }catch(e){showError(e.message)}finally{runBtn.querySelector('span').textContent='Run landscape analysis';updateBoundary()}
});
async function captureLeafletBasemap(){
  if(typeof html2canvas!=='function') throw new Error('Map screenshot library did not load. Refresh the page and try again.');
  // Capture the actual Esri tiles currently visible in Leaflet; omit draw handles and controls.
  const mapElement=map.getContainer();
  const canvas=await html2canvas(mapElement,{scale:2,useCORS:true,allowTaint:false,backgroundColor:'#101816',logging:false,ignoreElements:el=>el.classList?.contains('leaflet-control')||el.classList?.contains('leaflet-overlay-pane')||el.classList?.contains('leaflet-marker-pane')});
  const sw=L.CRS.EPSG3857.project(map.getBounds().getSouthWest());
  const ne=L.CRS.EPSG3857.project(map.getBounds().getNorthEast());
  // Retain exact rendered basemap pixels; polygon masking and flags are applied later in map coordinates.
  return {dataUrl:canvas.toDataURL('image/png'),bounds3857:[sw.x,sw.y,ne.x,ne.y]};
}
function showError(s){const el=document.getElementById('inlineError');el.textContent=s;el.hidden=false}
async function poll(){if(!currentJob)return;try{const r=await fetch(`/api/jobs/${currentJob}`);const j=await r.json();document.getElementById('progressBar').style.width=`${j.progress||0}%`;document.getElementById('progressText').textContent=`${j.progress||0}%`;document.getElementById('stageText').textContent=j.stage||j.status;document.getElementById('jobMessage').textContent=j.message||'Working…';document.getElementById('statusBadge').textContent=(j.status||'queued').toUpperCase();document.getElementById('logText').textContent=j.log||'Waiting for the pipeline…';
if(j.status==='completed'||j.status==='failed'){document.getElementById('referenceSection').hidden = j.status!=='completed';clearInterval(pollHandle);pollHandle=null;document.getElementById('jobTitle').textContent=j.status==='completed'?'Your landscape, decoded.':'Analysis needs attention';document.getElementById('jobMessage').textContent=j.status==='completed'?'Your study-area analysis and reports are ready.':(j.message||'The pipeline failed. Check the log below.');if(j.status==='completed')loadAutomaticReference();}
renderArtifacts(j.artifacts||[],j.status); }catch(e){document.getElementById('jobMessage').textContent='Cannot refresh status. The server may have restarted.'}}
function renderArtifacts(items,status){
  // Deduplicate on normalized relative path, then filename. This protects the UI even if an API returns duplicates.
  const seenPaths=new Set(), seenNames=new Set();
  const unique=items.filter(a=>{const path=String(a.path||a.name||'').replaceAll('\\','/').toLowerCase();const name=String(a.name||'').toLowerCase();if(seenPaths.has(path)||seenNames.has(name))return false;seenPaths.add(path);seenNames.add(name);return true;});
  const useful=unique.filter(a=>/interpretation_report\.md|summary\.json|method_comparison\.csv|method_comparison\.png|confidence\.tif|driver\.tif|missed\.tif|data_coverage\.png|layer_diagnostics\.csv|preprocessing_qa\.json|layer_trend_slopes\.tif|layer[1-4]_.*\.png|layer[1-4]_.*\.tif|flagged_region\.png|reference\.tif|trends_earth_reference\.tif/i.test(a.name));
  const signature=useful.map(a=>`${a.name}:${a.url}`).join('|')+`:${status}`;
  if(signature===lastArtifactSignature)return;lastArtifactSignature=signature;
  document.getElementById('outputCount').textContent=`${useful.length} FILE${useful.length===1?'':'S'}`;
  if(!useful.length){if(status==='failed')document.getElementById('resultsArea').innerHTML='<div class="empty-metrics">This run did not finish. Open the processing log to inspect the error.</div>';return}
  const area=document.getElementById('resultsArea');
  area.innerHTML='<div class="results-grid">'+useful.map(a=>`<article class="result-card"><div class="result-topline"><span class="result-type">${fileType(a.name)}</span><span class="result-icon">${fileIcon(a.name)}</span></div><h3>${escapeHtml(label(a.name))}</h3><p>${escapeHtml(a.name)}</p><a href="${escapeHtml(a.url)}">Download file <span>↗</span></a></article>`).join('')+'</div>';
  const csv=useful.find(a=>a.name.toLowerCase()==='method_comparison.csv');
  if(csv)loadMetrics(csv.url).catch(()=>showMetricFallback());
}
async function loadMetrics(url){
  const response=await fetch(url);if(!response.ok)throw new Error('Metrics file is not ready');const text=await response.text();
  const rows=parseCSV(text);if(rows.length<2)throw new Error('No comparison rows');
  const headers=rows[0].map(h=>h.trim());const records=rows.slice(1).map(row=>Object.fromEntries(headers.map((h,i)=>[h,(row[i]||'').trim()])));
  const fusion=records.find(r=>/four.?layer.?fusion/i.test(r.method||''))||records[records.length-1];
  const area=toNumber(fusion.area_km2),pixels=toNumber(fusion.pixels_flagged),f1=toNumber(fusion.F1),precision=toNumber(fusion.precision),recall=toNumber(fusion.recall),iou=toNumber(fusion.IoU);
  const methodCount=records.length;
  const metrics=[
    ['AREA FLAGGED',formatNumber(area,2,' km²'),`${fusion.method||'Fusion method'} · mapped footprint`],
    ['PIXELS FLAGGED',formatNumber(pixels,0),`${formatNumber(methodCount,0)} approaches compared`],
    ['F1 SCORE',isFinite(f1)?f1.toFixed(3):'N/A',isFinite(f1)?'Compared with reference labels':'Independent reference labels not supplied'],
    ['PRECISION / RECALL',isFinite(precision)&&isFinite(recall)?`${precision.toFixed(2)} / ${recall.toFixed(2)}`:'N/A',isFinite(precision)&&isFinite(recall)?'Reference-based evaluation':'Requires labelled reference raster']
  ];
  document.getElementById('metricsGrid').innerHTML=metrics.map(m=>`<article class="metric-card"><div class="metric-label">${escapeHtml(m[0])}</div><div class="metric-value">${escapeHtml(m[1])}</div><div class="metric-note">${escapeHtml(m[2])}</div></article>`).join('');
  document.getElementById('metricsSection').hidden=false;
}
function showMetricFallback(){document.getElementById('metricsGrid').innerHTML='<div class="empty-metrics">The method comparison file is present, but its metrics could not be read. Download the CSV to inspect the available values.</div>';document.getElementById('metricsSection').hidden=false}
function parseCSV(text){const rows=[];let row=[],field='',quoted=false;for(let i=0;i<text.length;i++){const c=text[i];if(c==='"'&&quoted&&text[i+1]==='"'){field+='"';i++}else if(c==='"'){quoted=!quoted}else if(c===','&&!quoted){row.push(field);field=''}else if((c==='\n'||c==='\r')&&!quoted){if(c==='\r'&&text[i+1]==='\n')i++;row.push(field);if(row.some(v=>v!==''))rows.push(row);row=[];field=''}else field+=c}if(field.length||row.length){row.push(field);rows.push(row)}return rows}
function toNumber(v){if(v===undefined||v===null||String(v).trim()==='')return NaN;const n=Number(v);return Number.isFinite(n)?n:NaN}
function formatNumber(v,digits,suffix=''){return Number.isFinite(v)?v.toLocaleString(undefined,{maximumFractionDigits:digits})+suffix:'N/A'}
function fileType(n){const ext=(n.split('.').pop()||'FILE').toUpperCase();return ext==='TIF'?'GEOSPATIAL RASTER':ext==='PNG'?'MAP / VISUAL':ext==='CSV'?'DATA TABLE':ext==='JSON'?'RUN METADATA':ext==='MD'?'FIELD REPORT':ext+' FILE'}
function fileIcon(n){if(/\.tif$/i.test(n))return '▦';if(/\.png$/i.test(n))return '▧';if(/\.csv$/i.test(n))return '≋';if(/\.json$/i.test(n))return '{ }';return '↗'}
function label(n){return ({'interpretation_report.md':'Interpretation report','summary.json':'Analysis summary','method_comparison.csv':'Method comparison table','method_comparison.png':'Method comparison chart','confidence.tif':'Degradation confidence raster','driver.tif':'Driver category raster','missed.tif':'Potential missed areas','data_coverage.png':'Data coverage chart','layer_diagnostics.csv':'Layer diagnostics','preprocessing_qa.json':'Preprocessing QA','layer1_ndvi.tif':'NDVI trend raster','layer2_recovery.tif':'Recovery-lag raster','layer3_restrend.tif':'RESTREND raster','layer4_sar.tif':'SAR trend raster','layer_trend_slopes.tif':'Trend slopes raster','layer1_ndvi_trend.png':'NDVI trend preview','layer2_recovery_lag.png':'Recovery-lag preview','layer3_restrend.png':'RESTREND preview','layer4_sar_structure.png':'SAR structure preview','l1_ndvi_trend_original_flagged.png':'L1 · Leaflet satellite image with NDVI flags','l2_recovery_lag_original_flagged.png':'L2 · Leaflet satellite image with recovery flags','l3_restrend_original_flagged.png':'L3 · Leaflet satellite image with RESTREND flags','l4_sar_structure_original_flagged.png':'L4 · Leaflet satellite image with SAR flags','flagged_region.png':'Synthesis comparison map','reference.tif':'Independent reference raster','trends_earth_reference.tif':'Trends.Earth regional reference raster'})[n]||n}
function escapeHtml(s){return String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}

// Display the automatically retrieved regional reference dataset and its provenance.
async function loadAutomaticReference(){
  const status=document.getElementById('automaticReferenceStatus');
  const links=document.getElementById('automaticReferenceLinks');
  status.textContent='Checking regional reference retrieval…'; links.innerHTML='';
  try{
    const response=await fetch(`/api/jobs/${currentJob}/reference-info`);
    const info=await response.json();
    if(info.evaluation_status==='metrics_calculated'){
      status.innerHTML=`<strong>Reference retrieved and evaluated.</strong><br>${escapeHtml(info.valid_coverage_percent??'—')}% of the analysis grid has reference coverage. Degraded pixels: ${formatNumber(Number(info.degraded_pixels),0)}; not-degraded pixels: ${formatNumber(Number(info.not_degraded_pixels),0)}. <br><small>${escapeHtml(info.reference_type||'External benchmark')}</small>`;
    }else if(info.evaluation_status==='reference_downloaded_metrics_unavailable'){
      status.innerHTML=`<strong>Reference raster retrieved, but binary metrics could not be calculated.</strong><br>${escapeHtml(info.evaluation_note||'The valid overlap may contain only one class.')}`;
    }else if(info.evaluation_status==='unavailable'){
      status.innerHTML=`<strong>Automatic reference data unavailable.</strong><br>${escapeHtml(info.error||info.note||'No further details were returned.')}<br><small>The four-layer screening outputs are still available, but are not reference-validated.</small>`;
    }else{
      status.textContent='No automatic reference information was recorded for this run.';
    }
    if(info.download_url){links.innerHTML=`<a href="${escapeHtml(info.download_url)}">Download aligned Trends.Earth reference raster ↗</a>`;}
    if(info.record_url){links.innerHTML+=` <a href="${escapeHtml(info.record_url)}" target="_blank" rel="noopener">Source dataset and documentation ↗</a>`;}
  }catch(error){status.textContent='Could not load reference metadata: '+error.message;}
}
