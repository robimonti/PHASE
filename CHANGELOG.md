# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Native Windows support for StaMPS stage (via [`pyccino/StaMPS`](https://github.com/pyccino/StaMPS) Windows fork), enabling end-to-end PSI processing on Windows 11.
- Automatic discovery of the StaMPS Windows installation via `StaMPS_CONFIG.ps1`, eliminating manual path configuration in `PHASE_StaMPS.mlapp`.
- Bidirectional Python interpreter sharing between PHASE and StaMPS: `setupPythonEnvironment` writes the chosen interpreter path to `%APPDATA%\PHASE\python.txt` on Windows so StaMPS can reuse it.
