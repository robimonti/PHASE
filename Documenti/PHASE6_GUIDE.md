# PHASE — archived detailed guide (PHASE 6 and development history)

This file preserves the former repository README for reference. For PHASE 7
installation and project organization, use the [current README](../README.md)
and [installer guide](../installer/README.md).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19605360.svg)](https://doi.org/10.5281/zenodo.19605360)
[![Latest release](https://img.shields.io/github/v/release/robimonti/PHASE)](https://github.com/robimonti/PHASE/releases/latest)

> [!IMPORTANT]
> PHASE 7 is being prepared for its first public release. Until a v7 release
> appears under [GitHub Releases](https://github.com/robimonti/PHASE/releases),
> the currently published v6.1.4 Windows installer remains the stable download.
> Do not use **Code → Download ZIP** as an installer.

> [!NOTE]
> The PHASE 7 installers described below are release candidates, not yet the
> published v6.1.4 installer.

**PHASE** (**P**ersistent scatterer **H**ighly **A**utomated **S**uite for **E**nvironmental monitoring) is a MATLAB-based software suite for automated InSAR Persistent Scatterer Interferometry (PSI) processing and advanced geospatial analysis. Built on the foundation of *snap2stamps* and *StaMPS*, PHASE introduces enhanced automation, user-friendly interactive map interfaces, and a powerful geospatial modeling module to interpret and visualize displacement time series, making it ideal for environmental and infrastructure monitoring.

PHASE 7 opens one hub with **Project**, **Preprocessing**, **StaMPS PSI**, and
**Displacement Modeling** sections. One installation manages projects stored
anywhere. The PHASE 6 applications remain in [`legacy`](legacy/README.md) for
provenance.

![Logo](https://github.com/user-attachments/assets/5bf0b784-c5e6-4e6c-8df5-2da8808263d3)

## SAR Satellites compatibility

- Sentinel-1 (from European Space Agency)
- COSMO-SkyMed (from Agenzia Spaziale Italiana - automatically supports both CSK and CSG generations)

## Required software

- SNAP (version 13.x is recommended)
- MATLAB (R2026a is recommended)
- StaMPS
- Python 3.10 or newer with `venv`; the installer prepares PHASE's Python packages

## Required OS
- *Windows 10/11*: guided installer and PSI runtime; tested on Windows.
- *macOS Apple Silicon*: guided DMG; preprocessing and StaMPS PSI tested with
  native Step 8 and TRAIN `a_linear`. Intel Macs are outside scope.
- *Linux x86_64/aarch64*: AppImage installer and native StaMPS/TRAIN build path
  prepared; an end-to-end PHASE 7 Linux run has **not** yet been tested.

## Installation and setup

> [!IMPORTANT]
Choose the installer for your OS from the matching
[GitHub Release](https://github.com/robimonti/PHASE/releases). The v7 assets
will be `install-phase.exe`, `PHASE-7-macos-arm64.dmg`, and a Linux AppImage.
MATLAB and SNAP are external prerequisites. StaMPS/TRAIN are installed or
built by the OS-specific setup. See the [installer guide](installer/README.md).

### Projects

Create projects from the hub. New projects use four numbered folders:
`01_INPUT` (optional user inputs), `02_PROCESSING_INTERNAL` (do not edit during
runs), `03_RESULTS` (final outputs), and `04_LOGS` (diagnostics). Project data
stays outside the PHASE installation. Existing projects keep their original
layout; PHASE does not rename paths that a processing session may depend on.
See the [project architecture](Documenti/PHASE_7_Architettura_Progetti.md).

Run `PHASE_Hub` in MATLAB to open the unified project workspace, or
`PHASE_Hub(projectRoot)` to open a specific PHASE 7 project. Preprocessing,
StaMPS and Model appear as sections of one window and load only when selected.
The standalone launchers remain available for diagnostics.

The OS-specific installer puts one engine outside project directories. The
hub's **Check for updates** button prepares a newer stable v7 engine from
GitHub; restarting from the installed shortcut applies it. The first
compatible update asset will ship with the first v7 release.

### Repository versus installed application

The GitHub **Code** page is the complete development repository. It includes
automated tests, release tooling, migration utilities and archived legacy files
so that scientific changes remain reproducible and maintainable. These are not
additional programs that a normal user needs to manage.

The current PHASE 7 installer sources remove development-only material
(`tests`, `legacy`, CI configuration and installer sources) from the installed
runtime. The installed hub launcher starts the editable MATLAB implementation
in the `engine` folder; projects live elsewhere. The published Windows 6.1.4
installer still presents three shortcuts; the v7 installer creates one.
The compiled installer is distributed as a GitHub **Release asset**, not
committed as a binary inside the source tree.

> [!NOTE]
> A detailed, step-by-step guide is available in the provided user manual. <br>
> Before using PHASE, please carefully read the entire manual!

### Installation

Use the installer attached to the release matching your version and OS. The
v7 Windows wizard creates one **PHASE** shortcut; the macOS DMG creates
`~/Applications/PHASE.app`; the Linux AppImage installs a PHASE launcher in the
user's application menu. [Installer details](installer/README.md).

### Historical manual setup (PHASE 6 only)

The following instructions document the older three-app PHASE 6 workflow.
They are **not** the installation procedure for PHASE 7; use the OS-specific
installer above instead.

1. **Install SNAP Software** <br>
   Download and install [SNAP 13.x](https://step.esa.int/main/download/snap-download/) from the European Space Agency website. <br> <br>
   Verify that the following mandatory SNAP plugin module is installed:
   - Sentinel-1 Toolbox <br>
   *(Note: Previous versions required multiple toolboxes like Optical or SMOS, but these are no longer needed).*<br>

   After installing SNAP, it is highly recommended to optimize your memory settings:
     - Edit `$HOME/snap/bin/gpt.vmoptions` and modify the `-Xmx` parameter according to your RAM (e.g., `-Xmx12G`).
     - Edit `$HOME/snap/etc/snap.properties` and add/verify:
          - `#snap.home=`
          - `#snap.userdir=`
          - `snap.jai.tileCacheSize = 1024`
          - `snap.jai.defaultTileSize = 512`

2. **Install Required Python Modules:** <br>
   Install [Python 3.x](https://www.python.org/downloads/) on your machine. Ensure Python is added to your system's PATH. The PHASE suite utilizes standard built-in Python libraries, so you only need to install the external Excel library. Run the following command in your terminal:
   ```bash
   pip install openpyxl requests asf_search shapely
   ```

3. **Install xterm (only Linux Users):** <br>
   Install xterm by running `sudo apt-get install xterm` in the terminal.

4. **Install StaMPS (only Linux Users):** <br>
   Install [StaMPS](https://homepages.see.leeds.ac.uk/~earahoo/stamps/) from the official GitHub repository.
   ```
   git clone https://github.com/dbekaert/StaMPS.git
   ```

5. **Install PHASE suite**
   - Download the latest release of the PHASE suite repository.
   - Move or extract the downloaded folder into your desired project directory.
   - Add the repository to the MATLAB path with `addpath(genpath(pwd))`.
   - Run `PHASE_Preprocessing` to open Module 1A immediately.
   - Tune the configurable parameters across the available tabs (including the interactive geographic map for AOI selection).
   - At the preprocessing handoff, PHASE opens `PHASE_StaMPS` for the generated `ASC_*` or `DSC_*` dataset. It may also be launched manually with `PHASE_StaMPS(datasetFolder)`.
   - Run `PHASE_Model` for geospatial analysis of the exported PS displacement time series.

The three launcher files are small production wrappers around the validated
editable engines. Internal package names retain the `_beta` suffix only for
backward compatibility with configurations and existing installations.

### Processing Steps

## Module 1: InSAR PSI Processing

1.	**Automated SAR Images Download:** <br>
Retrieve Sentinel-1 images via the integrated module through the Alaska SAR Facility APIs (thanks to Magnus and Johny). For COSMO-SkyMed, use the **Images** tab in the Cosmo-SkyMed panel to import your `.h5` files (they are copied into the `slaves` directory automatically).
2.	**Interactive AOI & Automated Master Selection:** <br>
Define your Area of Interest (AOI) by drawing a polygon directly on the integrated map or by importing supported geometry. Let PHASE automatically query the Open-Meteo historical weather API to select the optimal, driest master image for your stack.
3.	**Master & Slave Pre-Processing:** <br>
Automated splitting, precise orbit correction, coregistration, and interferogram formation. For Sentinel-1, optimal swaths and bursts are dynamically calculated from your AOI. Includes StaMPS export, average scene intensity computation, and local incidence angle/coherence calculations.
4.	**StaMPS Processing:** <br>
Automated data preparation, parameter definition, metadata auto-detection, and StaMPS PS analysis. Includes integration with TRAIN for GACOS tropospheric corrections and displacement time-series export. The optional **Export radar-grid diagnostics** control in the Export panel creates satellite-map figures and a candidate-level CSV documenting amplitude dispersion and the last StaMPS selection/weeding stage reached. The displayed grid is the SLC/interferometric sampling grid, not the physical SAR resolution.

## Module 2: Geospatial PSI Data Analysis

The geospatial module enhances PHASE by providing advanced interpolation and modeling of PS displacement time series, tailored for Sentinel-1 data but compatible with any SAR data in the same table format. It offers flexible processing options for environmental and infrastructure monitoring, with user-configurable parameters for experts and automated settings for beginners. Key features include:

- **Area of Interest (AOI) Definition**: Specify the AOI via shapefiles (preferred for arbitrary shapes) or bounding boxes, with automatic detection and transformation of geographic (WGS84) or projected (UTM) coordinates. PS outside the AOI are filtered to focus analysis.
- **Geometry Options**: Supports 1D modeling for linear features (e.g., roads, railways) using centerline interpolation and 2D modeling for expansive areas (e.g., volcanoes, large infrastructures) using grid-based meshes, with user-defined resolutions.
- **Temporal Modeling**: Independently interpolates each time series using cubic splines for outlier removal, trend, and periodic component modeling, followed by least-squares collocation for residuals. Optional nearest-neighbor interpolation extends results spatially for visualization.
- **Deterministic Spatio-Temporal Modeling**: Creates a continuous displacement field using temporal cubic splines for data cleaning and multi-dimensional splines for spatial interpolation, with per-epoch GIF visualizations.
- **Stochastic Spatio-Temporal Modeling**: Combines deterministic cubic spline-based cleaning with Least Squares Collocation for deformation modeling, incorporating covariance modeling for robust uncertainty estimates.
- **Alerts Analysis**: Evaluates PS stability via a thresholding procedure based on average velocity and cumulative displacement, assigning alert levels using standardized statistical criteria and exponential/linear alert scaling.

The module generates a comprehensive Excel report across multiple sheets:

- **General Sheet**: Summarizes project metadata (title, date, location, modeling approach), unwrapping parameters, PS statistics, and visualizations (logo, AOI map, PS plot).
- **Raw Displacement Sheet**: Presents adjusted raw displacement time series with coordinates and a figure of average scene displacement with velocity annotation.
- **Modeled Displacement Sheet**: Details modeled displacements with average velocity and a geoscatter plot of velocity across the AOI.
- **Uncertainty Sheet**: Provides uncertainty estimates with a geoscatter plot of time-averaged uncertainty.
- **Alerts Sheet**: Reports stability analysis with velocity, cumulative displacement, and global alert percentages.
- **Interpolation Sheet (Optional)**: Includes extrapolated time series at user-specified coordinates with individual displacement plots.

A shapefile and .mat file are generated in WGS84 (EPSG:4326) coordinates, including PS data, velocities, uncertainties, and risks for GIS compatibility.

## Possible Errors and Solutions
The procedure has been tested on SNAP 9.x and 13.x, Python 3.11, Python 3.13, Ubuntu 20.04, Windows 10, macOS Sequoia (15.1), and MATLAB 2025a/2026a. <br>
> [!TIP]
> Refer to the manual for solutions to common errors encountered during the StaMPS processing.

## SNAP version selection (Sentinel-1A/1B/1C/1D)

PHASE preprocessing supports both **SNAP 9.x** (legacy) and **SNAP 13.x** (recommended). Choose based on which Sentinel-1 satellites your dataset includes:

- **SNAP 9.x**: supports Sentinel-1A and Sentinel-1B only. Sentinel-1C (launched Dec 2024) and Sentinel-1D are **not** supported by the SNAP 9 product readers.
- **SNAP 13.x**: supports the full Sentinel-1A/1B/1C/1D constellation natively. **Required for any dataset containing S1C or S1D acquisitions.**

Point `GPTBIN_PATH` in your `project.conf` to the SNAP version you want PHASE to use:
```
GPTBIN_PATH = C:/Program Files/snap13/bin/gpt.exe   # SNAP 13 (S1A/B/C/D)
# GPTBIN_PATH = C:/Program Files/snap/bin/gpt.exe   # SNAP 9  (S1A/B only)
```

### Important: do not mix SNAP-9 and SNAP-13 .dim products

SNAP 13's `StampsExport` operator raises `NullPointerException` on tie-point grids written by SNAP 9 (the BEAM-DIMAP TPG format changed between majors). The error is silent — StaMPS later hangs mid-PSI without a clear diagnostic.

PHASE's `SEN_stamps_export.py` automatically detects this mismatch and prints an actionable warning before launching `gpt`. The fix is to **re-run the full preprocessing pipeline** (split + coregistration + interferogram) from the original `.SAFE.zip` files using the same SNAP version that will run `StampsExport`.

### SRTM 3Sec auto-cache

`StampsExport` has SRTM 3Sec hardcoded for the lat/lon geocoding output, regardless of the DEM chosen for coregistration. When the auto-download fails (offline, mirror down, or — on SNAP 13 — silently for some inputs), the geo files are partially corrupted and StaMPS hangs mid-PSI. PHASE pre-caches the required SRTM tiles into `%USERPROFILE%/.snap/auxdata/dem/SRTM 3Sec/` automatically before each `StampsExport` run, eliminating that failure mode.

This requires the project AOI to be present in `project.conf` as `LATMIN`, `LATMAX`, `LONMIN`, `LONMAX`. If the keys are absent the pre-cache is skipped (with a warning) and `StampsExport` falls back to the standard SNAP auto-download path.

## Verifying TRAIN on Windows

After the TRAIN Windows port, verify your install with these three checks.

### 1. Degradation path (TRAIN missing)

1. Open MATLAB. Run `which('aps_linear')`. Expected: empty string.
2. Launch `PHASE_StaMPS`. Tick "TRAIN atmospheric correction". Press Save, then Start.
3. Expected:
   - Warning id `StaMPS:phase:trainNotAvailable` in diary / `smoketest.log`.
   - The TRAIN checkbox STAYS TICKED (intentional — preserves intent for re-run).
   - Processing continues through STEP 1 → STEP 2 → export.
   - Output contains no `Atmosphere_*` columns.
   - Exit code 0.

### 2. Linear correction (`a_linear`)

1. Install TRAIN (Windows-patched fork): `git clone https://github.com/pyccino/TRAIN.git C:/TRAIN`.
2. In MATLAB: `addpath(genpath('C:/TRAIN/matlab')); savepath`.
3. Verify: `which('aps_linear')` returns `C:\TRAIN\matlab\aps_linear.m`.
4. Launch `PHASE_StaMPS`. Tick TRAIN. Set `tropo_method='a_linear'`. Save, Start.
5. Expected:
   - No degradation warning.
   - `aps_linear` runs (console output contains "loading the data").
   - Output contains `Atmosphere_a_linear_AOI_PS.mat` and `Atmosphere_a_linear_*.csv`.
   - Velocity values differ from a run with TRAIN unchecked.

> **Note on the Windows fork.** `pyccino/TRAIN` (default branch `main`) is forked from `dbekaert/TRAIN` at the audited commit `6c93feb` plus the following Windows-specific additions:
> - `get_gmt_version.m`: actionable error on Windows when GMT is not on PATH (the upstream loop manipulates Linux-only library env vars).
> - `aps_gacos_files.m`: replaces Unix `&` background launch with synchronous `system()` call on Windows (cmd.exe parses `&` differently).
> - `gacosDownloadDialog.m`: helper that shows the GACOS request parameters in a copy-paste dialog (called from PHASE StaMPS).
>
> Unix/Mac behavior is unchanged. Use upstream `dbekaert/TRAIN` directly on Linux/macOS if preferred.

### 3. GACOS correction (`a_gacos`) — optional, requires gacos.net data request

1. Same TRAIN install as above. Additionally install [GMT for Windows](https://www.generic-mapping-tools.org/download/) and ensure `C:\Program Files\GMT\bin` (or your install dir) is on PATH; verify with `gmt --version` in a fresh terminal.
2. Launch `PHASE_StaMPS`. Tick TRAIN. Set `tropo_method='a_gacos'`. Save, Start.
3. Expected:
   - A "Download GACOS maps" window opens (the gacos.net site and the `GACOS/` folder open automatically). It shows the request parameters (UTC, bounding box, dates) ready to copy into the form at gacos.net — select **Binary grid** as the file type.
   - Download the `.tar.gz` files from gacos.net, place them in `GACOS/` (do not extract), then press **Continue** in the window.
   - MATLAB extracts/distributes `.ztd` files.
   - Output contains `Atmosphere_a_gacos_AOI_PS.mat` and `Atmosphere_a_gacos_*.csv`.

## Updates
- *September 2026 — PHASE 6.1*: Added optional StaMPS radar-grid diagnostics: satellite-map figures of amplitude dispersion and candidate selection/weeding evolution, candidate-level CSV provenance, original SNAP lon/lat-grid geolocation, and multi-patch-safe reporting.
- *August 2026 — PHASE 6.0*: Promoted the standalone editable applications for preprocessing, StaMPS and geospatial modelling to production. Added the unified modern interface, integrated satellite maps and AOI drawing/import, in-app download and run monitoring, configurable processing controls, robust Windows runtime discovery, and the new installer layout. Archived the former `.mlapp` applications under `legacy`. Windows users should install from the release executable rather than cloning the repository.
- *June 2026*: Added the integrated download module for Sentinel-1. Completed the StaMPS porting to Windows; improved the StaMPS data export; created an installer for PHASE on Windows. Introduced the possibility to update the stack with newly available products abd update the pre-processing without re-starting from zero.
- *April 2026*: Added interactive geographic map GUI for automatic AOI sub-setting. Introduced meteorologically-aware master image selection using Open-Meteo API. Automated parameter metadata detection for StaMPS. Dropped legacy Python 2.7 support.
- *March 2026*: Introduced Module 2 for geospatial PSI data analysis with deterministic and stochastic modeling.
- *September 2024*: Added *macOS* compatibility to the preprocessing application and improved master error handling.
<img width="2100" height="1181" alt="GitHubUpdates" src="https://github.com/user-attachments/assets/f1fc37a6-67be-4770-a7d5-2a9b919ca4a3" />


## Planned updates
- Improve border constraints based on user selection.
- Introduce the handling of jumps in the displacement models.
- Add support for additional constellations.

## Acknowledgments
Special thanks to Jose Manuel Delgado Blasco and Dr. Michael Foumelis for the snap2stamps[^1] tool, and Prof. Andy Hooper for the StaMPS[^2] development. <br>

When using this software, please refer to:<br>
Monti, R., & Rossi, L. (2025). PHASE: a Matlab-based software for the DInSAR PS processing. Geodesy and Cartography, 51(2), 88–99. https://doi.org/10.3846/gac.2025.21995

[^1]: Foumelis, M., Delgado Blasco, J. M., Desnos, Y. L., Engdahl, M., Fernández, D., Veci, L. Lu, J. and Wong,
C. “SNAP - StaMPS Integrated processing for Sentinel-1 Persistent Scatterer Interferometry”. In
Geoscience and Remote Sensing Symposium (IGARSS), 2018 IEEE International, IEEE. <br>
[^2]: Hooper, A., A multi-temporal InSAR method incorporating both persistent scatterer and small baseline approaches, Geophys. Res. Lett., 35, L16,302, doi:10.1029/2008GL03465, 2008.
