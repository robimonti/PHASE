# PHASE Preprocessing 6.0 — implementation notes

`PHASE_Preprocessing.m` is the production entry point. It wraps the validated
`PHASE_Preprocessing_beta.m` implementation, whose internal name is retained
for backward compatibility. At runtime it neither loads nor copies
`PHASE_Preprocessing.mlapp` or `PHASE_StaMPS.mlapp`.

The formerly embedded class is generated as readable MATLAB source in
`PHASE_Preprocessing/+phase_preprocessing_beta/LegacyEngine.m`. Regenerate it
after stable MLAPP changes with:

```text
python tools/extract_preprocessing_beta_engine.py
```

The modern HTML interface covers the processing configuration, save/load,
execution, stop requests, live logs and data actions. Images is the first
workflow section: import is available for both constellations, while ASF
download and processed-stack update controls are shown only for Sentinel-1.

The Area of interest section contains its own interactive map. It displays only
the footprint read from the first valid Sentinel-1 ZIP or COSMO-SkyMed HDF5
product in `PHASE_Preprocessing/slaves`. The AOI can be drawn and edited beside
the four geographic fields; its bounding box updates the preprocessing
coordinates. Esri World Imagery supplies the satellite background and the
separate Boundaries and Places layer supplies city and place labels. Because
MATLAB `uihtml` does not support remote URL assets, the Python backend downloads
only the visible tiles into `phase_preprocessing_beta_ui/map_tiles`; the UI then
loads them as local supporting files. Cached tiles are reused. PHASE's vector
coastline remains as an offline fallback when map tiles are unavailable.

The first section now contains the complete image workflow without opening the
legacy MLAPP window:

- direct ZIP/HDF5 import into `PHASE_Preprocessing/slaves`;
- a live inventory table for the selected constellation;
- an independent download-selection map;
- the complete ASF filters, Earthdata login, result footprints, compatibility
  selection and download action;
- the existing-stack update search and download table.

ASF download and stack-update controls are present only when Sentinel-1 is the
selected constellation. Expanding the ASF filters or loading search results
extends the map column to the same height, so the selected footprints remain
visible beside the result table. Search results can be sorted by acquisition
date, path or frame; update results provide the same ordering choices.

Initial downloads and stack updates use the same background transfer manager.
It reports the overall percentage, current image/total, current-file byte
progress and filename without blocking MATLAB. Completed ZIPs are preserved,
unrelated files in `slaves` are never removed, interrupted `.part` files resume
on the next matching download, and transient network failures are retried up to
three times. **Force stop download** terminates the active transfer and leaves
the partial file available for resumption.

All of these controls call the existing PHASE search and authentication
backends; the extracted legacy engine remains hidden and is used only to run the
proven preprocessing callbacks.

## Processing configuration refinements

- Setup and Resources are one section. The Python/SNAP paths and SNAP CPU/cache
  settings are visible together.
- All processing controls are always visible; there is no advanced-controls
  visibility switch.
- When automatic master selection is enabled, the manual master date is
  disabled in the interface and remains ignored by validation.
- The first slave-processing step is shown as `number — description`. The
  sidebar range starts at that selected step and always terminates at step 6.
- Selecting a first step greater than 1 automatically disables master
  processing, because a resumed slave pipeline must reuse the existing master
  products. The Master and Slave sections show an explicit warning, and the
  MATLAB layer enforces the same rule for Sentinel-1 and COSMO-SkyMed.
- Interferogram/coregistration DEM selectors share a row with their external
  GeoTIFF paths. External paths are enabled only when `External DEM` is selected.
  The Sentinel-1 DEM-resampling selector occupies one half-row; cleanup controls
  begin on the following row.
- The historical `slaves_removal` option deletes the downloaded source
  ZIP/HDF5 files by removing and recreating `PHASE_Preprocessing/slaves`; it
  does not clean SNAP intermediates. The beta labels this explicitly, defers
  that source cleanup until the engine has returned successfully, and adds
  independent post-success cleanup switches for split/subset, coreg and ifg
  folders. No selected cleanup is performed after a failed run.

The Output products section can estimate a projected CRS from the AOI. It uses
the local WGS 84 UTM zone between 80°S and 84°N, including the Norway/Svalbard
zone exceptions, and polar stereographic CRSs outside that range. Existing MAT
configurations without the new `auto_epsg` field retain their manually selected
EPSG code.

Coherence and local incidence angle now have independent output switches. When
either is requested, the stable terrain-correction step runs; the unrequested
branch is not reprojected and is removed from the final output.

`+phase_preprocessing_beta/ProcessingEngine.m` subclasses the reproducibly
extracted legacy engine and corrects its GeoTIFF reprojection. The beta now
builds a target projected grid, inverse-projects its cell centres, resamples the
source raster in bounded memory blocks and atomically replaces the original
file. This is a true raster warp; it does not merely attach projected limits to
the unchanged source matrix.

## Live processing monitor

The beta no longer launches the generated master/slave scripts in separate CMD,
`xterm` or Terminal windows. It executes the same BAT/SH files through a silent
process bridge and streams merged stdout/stderr into Run monitor. SNAP
`gpt.exe` is also launched with the Windows no-console flag and its output is
forwarded while it runs rather than printed only at the end.

Run monitor remains fully visible during processing and includes the current
operation, milestone percentage, elapsed time and an ETA once enough progress
has been observed. The progress percentage is intentionally based on pipeline
milestones, not an invented per-pixel completion value. **Force stop now**
immediately terminates the active Windows CMD/Python/SNAP process tree. The
product being written can therefore be incomplete and must be recreated; PHASE
skips every optional cleanup after a forced stop.

On Windows a configured `...\snap\bin\gpt` value is resolved to the existing
`gpt.exe`; selecting the `bin` folder directly is also supported. When the
installer points to `pythonw.exe`, the beta uses the sibling `python.exe`
silently so Python and SNAP messages remain capturable in Run monitor.

After a successful export the beta explicitly creates the generated orbit/date
dataset directory. If `input_StaMPS.mat` is available, it is copied and
`PHASE_StaMPS_beta` opens normally. If it is not configured yet, preprocessing
still completes successfully, retains the new dataset directory and opens
`PHASE_StaMPS_beta` with an initial configuration. The StaMPS beta detects the
PHASE project and acquisition metadata when available; the user selects the
StaMPS installation folder and presses Save to create `input_StaMPS.mat`.

## Runtime independence

The two betas are self-contained with their editable backends:

- preprocessing uses `+phase_preprocessing_beta/LegacyEngine.m`;
- StaMPS uses `+phase_stamps_beta/runProcessing.m`;
- after preprocessing, the beta calls `PHASE_StaMPS_beta` directly and does not
  copy a legacy app into the generated `ASC_*`/`DSC_*` folder;
- the handoff invokes the canonical launcher from `PHASE_Preprocessing`, so an
  older manually copied launcher inside a dataset folder cannot shadow it.

The stable MLAPPs remain in the repository for production users and as the
provenance used by the developer extraction/regression tools. Those tools may
read the MLAPPs when regenerating or comparing source, but neither beta reads
them during normal execution.

On Windows, the preprocessing beta resolves Python 3 in this order: an
explicit executable selected in Setup, the interpreter recorded by the
installer in `%APPDATA%\PHASE\python.txt`, the Windows `py -3` launcher, and
finally the `python3`/`python` commands. Every candidate is version-checked
before it can run the map, ASF or preprocessing backends, so a legacy Python 2
installation earlier on `PATH` cannot be selected accidentally.

MATLAB-side configuration smoke test:

```matlab
phase_preprocessing_beta.selfTest
```

## Windows acceptance test

Test on a copy of the PHASE installation and of a small known dataset:

```matlab
cd('C:\path\to\PHASE')
addpath(genpath(pwd))
phase_preprocessing_beta.selfTest
phase_stamps_beta.selfTest
PHASE_Preprocessing_beta
```

1. Confirm that the satellite background and city labels appear in both the
   AOI and ASF download maps. On first view this can take a few seconds; reopen
   the section and confirm the cached background appears immediately.
2. Import one Sentinel-1 ZIP (or one COSMO-SkyMed HDF5), refresh the inventory,
   and confirm its footprint appears only after a valid file exists in
   `PHASE_Preprocessing/slaves`.
3. Draw and edit an AOI, save, close and reopen the beta, and confirm that the
   same bounding-box values are restored. Confirm Output products proposes the
   expected local UTM EPSG and that disabling automatic CRS enables manual EPSG.
4. Exercise ASF search/login/download only in an isolated test copy. Expand all
   filter groups and confirm the map remains visible beside the results; verify
   date/path/frame sorting.
5. Download two or three small test acquisitions and confirm that percentage,
   current image/total, bytes and the inventory advance during the transfer.
   Force-stop one download, restart the same selection and confirm that the
   `.part` file resumes. Repeat the progress/stop check for Stack update.
6. In Master processing, toggle automatic selection and confirm the manual date
   enables/disables accordingly. In Slave processing, verify the descriptive
   step range, paired DEM controls and each cleanup option. Select a first step
   greater than 1 and confirm master processing is disabled and the warning
   appears. Confirm all controls are immediately visible and CPU/cache are in
   Setup & resources.
7. Run a minimal known preprocessing case with every cleanup switch disabled.
   Confirm that no CMD window opens, master/Python/SNAP messages appear
   incrementally in Run monitor, and progress/elapsed time update during the
   run. Inspect the coherence/LIA/intensity GeoTIFF CRS and spatial alignment.
   Repeat with only one of coherence or LIA enabled.
8. Enable one intermediate cleanup option at a time on a disposable copy and
   confirm only the selected folder is removed after successful completion.
   On a separate disposable run, use **Force stop now**, verify that the active
   CMD/Python/SNAP tree ends immediately and discard the interrupted product.
9. Confirm preprocessing creates the generated orbit folder. If
   `input_StaMPS.mat` exists, confirm it opens `PHASE_StaMPS_beta`; otherwise
   confirm the StaMPS beta opens with detected initial values. Select the StaMPS
   installation folder, save to create the MAT file, then run the known StaMPS
   case first without TRAIN/GACOS and then with the atmospheric correction.
10. For the strongest independence check, in the copied installation only,
   temporarily rename both legacy `.mlapp` files to `.mlapp.disabled` and
   repeat the launch and smoke tests. Restore their names afterwards.
