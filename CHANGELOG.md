# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [7.1.8] - 2026-11-04
- Consolidate the validated custom-model workflow: exclusive expression/function selection, shapefile map preview, clearer coordinate/preview controls, covariance ordering and output sections.
- Replace hub alerts wording with relative anomalies; remove temporary UI actions while preserving the validated scientific calculations.
- Align application and all installers to 7.1.8.

## [7.1.7] - 2026-11-03
- Remove temporary example UI; choose whether to restore the automatic custom residual field, with target-consistent uncertainty and component exports.
- Enable native temporal checks by default, including local PS trends and identifiable time-space interactions; preserve atmospheric-offset semantics.
- Fix degree-three 1D STC observation reconstruction; align app and all installers to 7.1.7.

## [7.1.6] - 2026-11-03
- Simplify custom-model controls and add an estimated temporal quadratic example. Preserve existing advanced profile settings and expose their active status.
- Align app and installers to 7.1.6.

## [7.1.5] - 2026-11-03
- Separate temporal trend selection from spatial prediction validation; refine regularisation among temporally supported mean fields.
- Add independent-noise, oscillatory-signal and valid long-range covariance regression checks.
- Align application and installer versions to 7.1.5.
- Prepare 7.1.4: group covariance controls, consolidate editable figures and automatically validate residual mean fields in custom modelling. Preserve native outlier/covariance estimators and record unsupported fits for review.

## [7.1.3] - 2026-11-03
- Restore native covariance figures and organize custom outputs in files/customModel.
- Add optional residual temporal diagnostics/correction and bounded covariance lengths across automatic/custom DET/STC. Preserve original outlier filtering.
- Propagate the additional custom temporal fit with shared-data cross covariance.
- Autofit numeric report columns, wrap sidebar paths and space the AOI source controls.
- Fix skipped-candidate preallocation in automatic DET2D; align app/installer versions.

## [7.1.2] - 2026-11-02
- Export custom STC empirical covariance/kernel plots as PNG and FIG, with
  contributing pair counts, marginal fits and the actual common-amplitude
  prediction kernel; export residual-mean diagnostic plots for custom DET/STC.
- Recover covariance plots from an existing STC JSON without refitting.
- Enter query-point IDs and WGS84 coordinates in the GUI or load TXT/CSV;
  GUI points are saved to the same TXT format and use the same processing path.
- Export dated custom query series with components and conditional uncertainty;
  fix X/Y query detection and normalise numeric station IDs for report exports.
- Preserve input coordinates/times and query predictions in custom MAT results.
- App and installer labels aligned to 7.1.2; automatic scientific engines and
  outlier rejection remain unchanged.

## [7.1.1] - 2026-11-02
- Simplified Automatic/Custom selection and moved the AOI source choice above
  the map; shapefile mode no longer requires opening the drawing map.
- Custom model time starts at the first acquisition, with no separate offset.
- STC memory budget is expressed in GiB with an explicit matrix-RAM estimate.
- Custom progress distinguishes outlier filtering, model fitting, covariance
  fitting/factorisation and grid/PS/query prediction with uncertainty.
- Equivalent Cholesky-whitened uncertainty equations and larger bounded
  prediction blocks remove repeated dense covariance products. Automatic
  scientific functions and the original outlier filter remain unchanged.
- Clarified that the output temporal step is in days, including predictions
  between acquisitions. Private student input is excluded from Git.
- App/installer development versions aligned to 7.1.1.
- Automatic spatial centering for the normal custom workflow; coordinate
  overrides are grouped in an optional section, expanded for non-defaults.

## [7.1.0] - 2026-10-30
- Optional custom deterministic models for spatio-temporal DET/STC in 1D/2D,
  with fixed or bounded estimated parameters and MATLAB function support.
- Separate custom entry points preserve the automatic scientific engines and
  copy their outlier rejection verbatim.
- Two-stage fitting, explicit shared-data covariance propagation and
  reproducible component/parameter exports. See
  [the custom-model guide](Documenti/PHASE_7_1_CUSTOM_MODELS.md).
- App and installer labels aligned to 7.1.0; custom models released for real-project testing.

### Added
- Unified PHASE 7 hub and one-installation/many-project workflow.
- Apple Silicon installer and tested native StaMPS/TRAIN PSI path.
- Linux AppImage installer with pinned local StaMPS/TRAIN build preparation.
- In-app v7 engine update flow and draft release packaging workflow.
- Sequential four-folder layout for new projects; existing layouts preserved.

### Fixed
- StaMPS Python module resolution on macOS and SNAP memory guards across OSes.
- Incorrect wrapped-phase completion text after an unwrapped export.

## [6.1.1] - 2026-09-28
### Fixed
- Prevented the Windows installer from hanging or leaving orphaned Git
  processes when installing to mapped or UNC network drives.
- Clone repositories locally with bounded subprocesses, then copy them to the
  chosen destination with retry-limited `robocopy`.

## [6.1.0] - 2026-09-17
### Added
- Optional StaMPS radar-grid diagnostic figures for initial amplitude
  dispersion and candidate selection/weeding evolution.
- Candidate-level CSV and JSON provenance export using the original SNAP
  longitude/latitude grids, with an explicit interpolation fallback.
- Multi-patch-safe stage reporting and merge-resampling metadata.

## [6.0.0] - 2026-08-01
### Added
- Native Windows support for StaMPS stage (via [`pyccino/StaMPS`](https://github.com/pyccino/StaMPS) Windows fork), enabling end-to-end PSI processing on Windows 11.
- Automatic discovery of the StaMPS Windows installation via `StaMPS_CONFIG.ps1`, eliminating manual path configuration in `PHASE_StaMPS.mlapp`.
- Bidirectional Python interpreter sharing between PHASE and StaMPS: `setupPythonEnvironment` writes the chosen interpreter path to `%APPDATA%\PHASE\python.txt` on Windows so StaMPS can reuse it.
