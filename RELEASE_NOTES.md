# PHASE v6.1.4

PHASE v6.1.4 fixes an upstream StaMPS Step 3 defect exposed by non-zero
`gamma_stdev_reject` and preserves two decimal places in the weeding standard
deviation field.

## StaMPS candidate-selection hotfix

- PHASE detects and repairs the known `sfprintf` typo in the locally installed
  `StaMPS/matlab/ps_select.m` before processing begins. The original is saved
  once as `ps_select.m.phase-original`.
- The failure affected only configurations with `gamma_stdev_reject > 0`; it
  was not caused by Sentinel-1 data, preprocessing output, or parameter type
  conversion.
- `weed_standard_dev` now supports and preserves up to two decimal places, so
  values such as `0.95` reach StaMPS unchanged.

## Previous v6.1.3 changes

PHASE v6.1.3 makes long Sentinel-1 ASF download queues recoverable and prevents
the preprocessing interface from remaining permanently busy after a failed
image or an optional inventory-refresh error.

## ASF downloader reliability

- A failed image no longer blocks the remaining queue. PHASE records the
  failure and continues with every subsequent selected acquisition.
- Each image is retried up to five times with exponential backoff and support
  for ASF `Retry-After` throttling responses.
- Earthdata authentication and the HTTPS session are refreshed between retries
  during long-running downloads.
- Completed ZIP files and resumable `.part` files are preserved.
- Failed or stopped transfers now expose a **Retry pending downloads** button;
  already completed files are detected and skipped automatically.
- Download manifest, progress, result and Python traceback are retained under
  `engine/downloadasf/last_download*` for diagnosis instead of being deleted.
- ZIP inventory and footprint-refresh failures can no longer interrupt transfer
  finalisation or leave the UI controls disabled.
- Queue progress reaches the end even when individual images fail, while the
  final result still reports every failed scene explicitly.

## Previous v6.1.2 changes

PHASE v6.1.2 fixes the compiled Windows installer's repository downloader.

## Installer hotfix

- Replaced the unreliable `Start-Process` wrapper with a .NET process runner
  that preserves the real exit code and captures Git/Robocopy diagnostics when
  compiled with PS2EXE.
- If anonymous `git clone` is rejected or interrupted by a proxy, credential
  manager, antivirus, or HTTP transport issue, the installer automatically
  downloads the identical public GitHub branch archive and continues.
- Repository acquisition still happens in local temporary storage before the
  completed files are copied to mapped, UNC, or SMB destinations.
- Timeouts still terminate the entire process tree, preventing orphaned Git
  processes.
- Existing archive-based installations can now be refreshed even though they
  intentionally do not contain a `.git` directory.

## Previous v6.1.1 changes

PHASE v6.1.1 fixes Windows installation on mapped and network drives and
includes the optional, traceable diagnostic export introduced in v6.1.

## Installer reliability fix

- Git repositories are cloned into a local temporary staging directory before
  being copied to the selected installation path. Git no longer performs
  checkout operations directly on mapped/UNC/SMB storage.
- Every Git and copy subprocess has a hard timeout. On timeout PHASE terminates
  the complete child process tree, including `git-remote-https.exe`, so closing
  the installer cannot leave orphaned Git processes behind.
- Repository downloads are shallow and non-interactive, reducing installation
  time and preventing hidden credential prompts.
- Completed installations can be refreshed from local staging without deleting
  untracked PHASE project data already present below the engine folder.
- Failed first-time copies remove their partial destination and local staging
  data before returning an actionable error.

PHASE v6.1 added an optional, traceable diagnostic export for the StaMPS
candidate-selection workflow while retaining the complete standalone PHASE 6
application suite introduced in v6.0.

## Windows installation — important

**Windows users should download and run `install-phase.exe` attached to this
release. Do not clone the repository for a normal Windows installation.**

The installer creates and configures the complete PHASE working environment,
including StaMPS, TRAIN, the mandatory native Windows executables, Python,
MATLAB paths, and shortcuts for all three applications. A repository clone is
intended only for PHASE development and does not perform this setup.

After installation, launch the applications using:

- `PHASE Preprocessing.lnk`
- `PHASE StaMPS.lnk`
- `PHASE Model.lnk`

The editable MATLAB source remains available in the visible `engine` folder.

## New in v6.1

- Added **Export radar-grid diagnostics** to panel 11, Export, in PHASE StaMPS.
- Added a satellite-map figure of initial candidates coloured by amplitude
  dispersion `D_A`.
- Added a second satellite-map figure showing whether each initial candidate
  was rejected by `ps_select`, rejected during weeding, or survived patch-local
  weeding.
- Added a candidate-level CSV containing patch provenance, zero-based
  azimuth/range indices, coordinates, `D_A`, the effective threshold and the
  maximum processing stage reached.
- Added JSON metadata documenting grid meaning, geolocation source, merge
  resampling, class counts and grid decimation.
- Uses the original full-resolution SNAP `.lon`/`.lat` grids already referenced
  by `psclonlat.in`, with automatic byte-order and candidate-centre validation.
- Supports multiple StaMPS patches without incorrectly equating patch-local
  weeding survivors to the merged root `ps2` product.
- Reads the effective candidate threshold from the run's `selpsc.in` rather
  than relying on a possibly stale interface value.
- Retains an explicitly reported local interpolation fallback when original
  SNAP geolocation rasters are unavailable.

The displayed radar grid is the SLC/interferometric sampling grid. It must not
be interpreted as the physical SAR resolution or point-spread function.

## PHASE 6 application suite

- Three standalone, editable MATLAB applications with no runtime dependency on
  the former App Designer `.mlapp` files.
- Unified clean interface for preprocessing, StaMPS and geospatial modelling.
- Integrated Sentinel-1 search, download, stack update and progress monitoring.
- Interactive satellite maps for footprint inspection and polygon AOI drawing.
- In-app run monitoring, progress, elapsed time, ETA and hard-stop controls.
- Windows-native StaMPS/TRAIN discovery and automatic installation of required
  executables, including `snaphu.exe` and `triangle.exe`.
- Configurable StaMPS step range, temporal windows and scientific parameters.
- Improved geospatial modelling, temporal-only reporting, shapefile AOIs,
  GeoSPLINTER execution and robust interpolation of duplicate sample points.
- Former `.mlapp` applications archived under `legacy` for provenance only.

## Other platforms

Linux and macOS users, and developers who need the source tree, can use the
manual installation instructions in the repository README.
