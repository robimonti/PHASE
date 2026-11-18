# Preliminary PSI suitability

After preprocessing, open **PSI suitability** in the Preprocessing sidebar.
The map is generated automatically before optional SNAP-product cleanup.
**Build quality map** also evaluates an existing, complete SNAP export without
running preprocessing again. It reads the saved processing configuration.

The same layers are available on the optional PSI sub-area selection map.
Select a layer, adjust opacity, and move over the map to inspect its values and
the reasons for a screening class. The map never changes the preparation mask,
PS selection, weeding, unwrapping or displacement modelling.

## Layers

| Layer | Meaning |
|---|---|
| Preliminary suitability | Favourable, Intermediate, Unfavourable or Not evaluated; heuristic screening classes, not PS probabilities. |
| Terrain mapping risk | Fraction of sampled ground positions in a cell with terrain folding, back-facing slopes or an inconsistent terrain/radar/coordinate round trip. |
| Amplitude stability | Fraction of sufficiently observed radar pixels with amplitude dispersion below 0.4 after acquisition-gain normalisation. |
| Stack coverage | Mean fraction of acquisitions with finite, nonzero complex samples at the observed radar pixels. |
| Interferogram coherence | Median, across retained interferograms, of each cell's spatial mean coherence. This is not StaMPS temporal coherence. |

The stability colour scale saturates at 5%; hover values and GeoTIFF values
retain the actual fraction. Low stable-pixel fractions do not mean that every
individual target is unsuitable. Coherence depends on the spatial window,
temporal and perpendicular baselines, and scene properties. The report records
the interferograms' temporal baselines.

The geographic grid is at least 50 m for CSK and 100 m for Sentinel-1, and may
be coarser for large AOIs. Pixel counts are radar samples, not ground-area
fractions. The layer resolution is an analysis scale, not a claimed position
accuracy or density of recoverable PS.

## Geometry availability and limits

Terrain screening uses SNAP's installed DEM/datum handling and orbit model.
Before accepting the model, PHASE projects exported coordinates and their
ellipsoidal heights back into radar coordinates. The 95th percentile residual
must be at most one pixel, and no residual may exceed three pixels. Missing or
incompatible metadata leaves the geometric indicator **not evaluated**.

The terrain grid is bounded for memory and runtime. Its actual sampling
interval is recorded. The screen detects local folding and slopes facing away
from the radar, and checks whether a ground-to-radar-to-exported-coordinate
path returns within the declared tolerance. It is **not a validated exhaustive
layover/shadow mask**, does not trace all cast shadows, and cannot certify
individual roof heights, absolute geolocation or successful PS recovery.
Favourable means favourable evidence at the analysis scale, not guaranteed
observability. These screening classes have not been calibrated as
sensor-independent probabilities.

## Transparent screening rules

At least five acquisitions and twenty geolocated radar samples per cell are
required. Amplitude dispersion requires observations in at least 80% of the
stack and at least five dates. Each acquisition is divided by its mean finite,
nonzero amplitude inside the AOI. Population standard deviation is used.

- Unfavourable evidence: coverage below 80%, stable-pixel fraction below 1%,
  median interferogram coherence below 0.2, or terrain-risk fraction at least 20%.
- Favourable evidence: coverage at least 95%, stable-pixel fraction at least 5%,
  coherence at least 0.4, and terrain-risk fraction below 5%, with every
  indicator available.
- Intermediate: sufficient observations and all indicators available, without
  satisfying either of the preceding rules.
- Missing evidence or insufficient observations never produces a favourable
  class. Poor measured evidence can still produce an unfavourable class while
  another indicator is unavailable; hover inspection reports both reasons.

These are explicit initial screening rules, not fitted likelihoods. Terrain
risk is never cancelled out by high amplitude stability or coherence.

## Export and reproducibility

Files are stored in `INSAR_<master>/quality/`:

- `suitability.tif`, `geometry.tif`, `stability.tif`, `coverage.tif`,
  `coherence.tif`: WGS84 GeoTIFFs, north-up, NaN for unavailable/outside data.
- `quality.json`: metric definitions, thresholds, normalisation factors,
  source-file identities, geometry validation and per-cell reasons. Class 0
  means Not evaluated; classes 1/2/3 mean Favourable/Intermediate/Unfavourable.
- `figures/PSI_quality.png` and `figures/PSI_quality.fig`: overview and editable
  MATLAB figure.
- `geometry.log`: diagnostic output when a geometry projection was attempted.

Source rasters are never changed. Statistics are streamed in bounded blocks.
The report becomes visible only after all layers have been written. Changes
to retained source data invalidate the cached map. Optional SNAP products
removed by the configured cleanup do not invalidate an already generated
report; a rebuild after cleanup may lack their coherence or terrain evidence.

Background references: [EGMS algorithm documentation](https://library.land.copernicus.eu/products/European_Ground_Motion_Service_Algorithm_Theoretical_Basis_Document_v4.html),
[ESA SAR Basics](https://step.esa.int/docs/tutorials/S1TBX%20SAR%20Basics%20Tutorial.pdf).
