# Land Degradation Detection — Interpretive Results

**Study period:** 2001–2025
**Valid analysis pixels:** 615
**Estimated valid area:** 38.44 km²
**Pixels flagged by at least one layer:** 172 (27.97% of valid pixels)

## How to interpret the four layers

- **L1 NDVI trend:** persistent decline in annual vegetation greenness. This is the baseline and is relatively easy to explain, but it can miss structural changes and resilience loss.
- **L2 recovery lag:** increasing time for vegetation anomalies to recover after rainfall-defined dry years. Treat as a resilience indicator; event count and missing observations matter.
- **L3 RESTREND-style residual:** declining NDVI after a simple linear rainfall adjustment. It suggests a decline not explained by this model's rainfall relationship, not proof of a human cause.
- **L4 SAR structure:** declining relative Sentinel-1 VH backscatter in the selected seasonal window. This can add structural evidence but can also respond to soil moisture and surface conditions.

## Comparison with alternative approaches

| Method | Flagged area (km²) | What it adds |
|---|---:|---|
| NDVI only | 6.00 | Simple greenness-trend baseline; can miss resilience, rainfall-adjusted and radar-only signals. |
| NDVI plus RESTREND | 6.06 | Adds a simple rainfall-adjusted residual trend to the NDVI baseline. |
| Optical plus resilience | 6.06 | Adds recovery-lag evidence but remains dependent on optical greenness. |
| Four layer fusion | 10.75 | Combines greenness, recovery, rainfall-adjusted residuals and SAR structural-change evidence. |

Four-layer area not flagged by NDVI-only baseline: **4.75 km²**.
This is a difference in mapped detection footprint, not an accuracy gain. A larger footprint can include true positives and false positives.

Independent reference labels were supplied, so classification metrics were calculated on aligned valid pixels:

| Method | Precision | Recall | F1 | IoU | Balanced accuracy | Reference pixels |
|---|---:|---:|---:|---:|---:|---:|
| NDVI only | 0.990 | 0.438 | 0.607 | 0.436 | 0.717 | 477 |
| NDVI plus RESTREND | 0.990 | 0.442 | 0.611 | 0.440 | 0.719 | 477 |
| Optical plus resilience | 0.990 | 0.442 | 0.611 | 0.440 | 0.719 | 477 |
| Four layer fusion | 0.842 | 0.613 | 0.709 | 0.550 | 0.758 | 477 |

Use these scores to compare alternatives on the same reference pixels. Check class balance and spatial/temporal independence before making a superiority claim.

## Confidence score

The confidence raster counts how many of the four layer flags are positive (0–4). It is an evidence-agreement score, not a probability. Layer outputs can be correlated because several depend on NDVI or rainfall.

## Limitations and recommended validation

1. Thresholds are configurable and should be sensitivity-tested; results may change with alpha, slope cut-offs, minimum event counts and the study period.
2. The recovery-lag method uses simplified drought and recovery definitions. Inspect representative time series and compare events with known drought records.
3. The RESTREND implementation is a simple rainfall-adjusted linear residual trend; temperature, land cover, irrigation and other drivers are not modelled. Annual rainfall is used only when all 12 monthly totals are present.
4. Sentinel-1 VH change is not uniquely attributable to canopy degradation; soil moisture, roughness and acquisition geometry can affect it.
5. If independent reference labels are unavailable, do not claim higher accuracy than alternatives. Use field data, expert-labelled samples or an independent degradation product for validation.
6. Compare alternatives on the same valid pixels and period. Report precision, recall, F1, IoU and area-adjusted accuracy only after reference labels exist.

## Output guide

- `flagged_region.png`: side-by-side layer flags, synthesis, heuristic driver, NDVI-only omissions and area comparison.
- `confidence.tif`: 0–4 layer-agreement count.
- `driver.tif`: heuristic category only; not causal attribution.
- `layer_trend_slopes.tif`: continuous slope/effect-size surfaces for all four layers.
- `layer_diagnostics.csv`: pixel counts, mapped area, valid slope support and median slope per layer.
- `method_comparison.csv`: comparable detection footprints and, when reference labels are configured, precision/recall/F1/IoU metrics.
- `summary.json`: machine-readable metrics and diagnostics.

Generated: 2026-11-02T13:03:10