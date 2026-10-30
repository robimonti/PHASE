# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
