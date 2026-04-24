# PHASE (Persistent scatterer Highly Automated Suite for Environmental monitoring)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19605360.svg)](https://doi.org/10.5281/zenodo.19605360)

**PHASE** (**P**ersistent scatterer **H**ighly **A**utomated **S**uite for **E**nvironmental monitoring) is a MATLAB-based software suite for automated InSAR Persistent Scatterer Interferometry (PSI) processing and advanced geospatial analysis. Built on the foundation of *snap2stamps* and *StaMPS*, PHASE introduces enhanced automation, user-friendly interactive map interfaces, and a powerful geospatial modeling module to interpret and visualize displacement time series, making it ideal for environmental and infrastructure monitoring.

![Logo](https://github.com/user-attachments/assets/5bf0b784-c5e6-4e6c-8df5-2da8808263d3)

## SAR Satellites compatibility

- Sentinel-1 (from European Space Agency)
- COSMO-SkyMed (from Agenzia Spaziale Italiana - automatically supports both CSK and CSG generations)

## Required Software
- SNAP (version 9.x is mandatory)
- MATLAB
- StaMPS
- Python 3.x with `openpyxl` (for geospatial module report generation)

## Required OS
- *Linux*: supported end-to-end, including StaMPS, for the complete PSI pipeline.
- *Windows*: supported end-to-end for the SNAP-based workflow (preprocessing,
  StaMPS, and geospatial analysis) via the Windows-native StaMPS fork
  [`pyccino/StaMPS`](https://github.com/pyccino/StaMPS). PHASE auto-discovers
  the StaMPS Windows install through `StaMPS_CONFIG.ps1` and shares the Python
  interpreter with StaMPS via `%APPDATA%\PHASE\python.txt`.
- *macOS*: supports preprocessing and geospatial analysis modules.

## Installation and Setup

> [!NOTE]
> A detailed, step-by-step guide is available in the provided user manual. <br>
> Before using PHASE, please carefully read the entire manual!

### Preliminary Steps
1. **Install SNAP Software** <br>
   Download and install [SNAP 9.x](https://step.esa.int/main/download/snap-download/) from the European Space Agency website. <br> <br>
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
   pip install openpyxl

3. **Install xterm (only Linux Users):** <br>
   Install xterm by running `sudo apt-get install xterm` in the terminal.

5. **Install StaMPS:** <br>
   Install [StaMPS](https://homepages.see.leeds.ac.uk/~earahoo/stamps/) from the official GitHub repository.
   ```
   git clone https://github.com/dbekaert/StaMPS/releases/tag/v4.1-beta
   ```

6. **Install PHASE suite**
   - Download the latest release of the PHASE suite repository.
   - Move or extract the downloaded folder into your desired project directory.
   - Execute the PHASE_Preprocessing.mlapp MATLAB application.
   - Tune the configurable parameters across the available tabs (including the interactive geographic map for AOI selection).
   - Once the preprocessing is complete, execute the PHASE_StaMPS.mlapp MATLAB application (on Linux, it will open automatically upon completion).
   - For geospatial analysis, run PHASE_model.mlapp to process the final displacement time series.

### Processing Steps

## Module 1: InSAR PSI Processing

1.	**Automated SAR Images Download:** <br>
Retrieve Sentinel-1 images via the generated Python script from the Alaska SAR Facility, or manually place your COSMO-SkyMed .h5 files into the slaves directory.
2.	**Interactive AOI & Automated Master Selection:** <br>
Define your precise Area of Interest (AOI) by drawing a bounding box directly on the GUI's geographic map interface. Let PHASE automatically query the Open-Meteo historical weather API to select the optimal, driest master image for your stack.
3.	**Master & Slave Pre-Processing:** <br>
Automated splitting, precise orbit correction, coregistration, and interferogram formation. For Sentinel-1, optimal swaths and bursts are dynamically calculated from your AOI. Includes StaMPS export, average scene intensity computation, and local incidence angle/coherence calculations.
4.	**StaMPS Processing:** <br>
Automated data preparation, parameter definition, metadata auto-detection, and StaMPS PS analysis. Includes integration with TRAIN for GACOS tropospheric corrections, exporting displacement time series in Excel format.

## Module 2: Geospatial PSI Data Analysis

The geospatial module enhances PHASE by providing advanced interpolation and modeling of PS displacement time series, tailored for Sentinel-1 data but compatible with any SAR data in the same table format. It offers flexible processing options for environmental and infrastructure monitoring, with user-configurable parameters for experts and automated settings for beginners. Key features include:

- **Area of Interest (AOI) Definition**: Specify the AOI via shapefiles (preferred for arbitrary shapes) or bounding boxes, with automatic detection and transformation of geographic (WGS84) or projected (UTM) coordinates. PS outside the AOI are filtered to focus analysis.
- **Geometry Options**: Supports 1D modeling for linear features (e.g., roads, railways) using centerline interpolation and 2D modeling for expansive areas (e.g., volcanoes, large infrastructures) using grid-based meshes, with user-defined resolutions.
- **Temporal Modeling**: Independently interpolates each time series using cubic splines for outlier removal, trend, and periodic component modeling, followed by least-squares collocation for residuals. Optional nearest-neighbor interpolation extends results spatially for visualization.
- **Deterministic Spatio-Temporal Modeling**: Creates a continuous displacement field using temporal cubic splines for data cleaning and multi-dimensional splines for spatial interpolation, with per-epoch GIF visualizations.
- **Stochastic Spatio-Temporal Modeling**: Combines deterministic cubic spline-based cleaning with Least Squares Collocation for deformation modeling, incorporating covariance modeling for robust uncertainty estimates.
- **Stability Analysis**: Evaluates PS stability via a thresholding procedure based on average velocity and cumulative displacement, assigning risk levels using standardized statistical criteria and exponential/linear risk scaling.

The module generates a comprehensive Excel report across multiple sheets:

- **General Sheet**: Summarizes project metadata (title, date, location, modeling approach), unwrapping parameters, PS statistics, and visualizations (logo, AOI map, PS plot).
- **Raw Displacement Sheet**: Presents adjusted raw displacement time series with coordinates and a figure of average scene displacement with velocity annotation.
- **Modeled Displacement Sheet**: Details modeled displacements with average velocity and a geoscatter plot of velocity across the AOI.
- **Uncertainty Sheet**: Provides uncertainty estimates with a geoscatter plot of time-averaged uncertainty.
- **Alerts Sheet**: Reports stability analysis with velocity, cumulative displacement, and global risk percentages.
- **Interpolation Sheet (Optional)**: Includes extrapolated time series at user-specified coordinates with individual displacement plots.

A shapefile and .mat file are generated in WGS84 (EPSG:4326) coordinates, including PS data, velocities, uncertainties, and risks for GIS compatibility.

## Possible Errors and Solutions
The procedure has been tested on SNAP 9.x, Python 2.7, Python 3.11, Ubuntu 20.04, Windows 10, macOS Sequoia (15.1), and MATLAB 2025a. <br>
> [!TIP]
> Refer to the manual for solutions to common errors encountered during the StaMPS processing.

## Verifying TRAIN on Windows

After the TRAIN Windows port, verify your install with these three checks.

### 1. Degradation path (TRAIN missing)

1. Open MATLAB. Run `which('aps_linear')`. Expected: empty string.
2. Launch `PHASE_StaMPS.mlapp`. Tick "TRAIN atmospheric correction". Press Save, then Start.
3. Expected:
   - Warning id `StaMPS:phase:trainNotAvailable` in diary / `smoketest.log`.
   - The TRAIN checkbox STAYS TICKED (intentional — preserves intent for re-run).
   - Processing continues through STEP 1 → STEP 2 → export.
   - Output contains no `Atmosphere_*` columns.
   - Exit code 0.

### 2. Linear correction (`a_linear`)

1. Install TRAIN at the audited commit: `git clone https://github.com/dbekaert/TRAIN.git C:/TRAIN && cd C:/TRAIN && git checkout 6c93feb`.
2. In MATLAB: `addpath(genpath('C:/TRAIN/matlab')); savepath`.
3. Verify: `which('aps_linear')` returns `C:\TRAIN\matlab\aps_linear.m`.
4. Launch `PHASE_StaMPS.mlapp`. Tick TRAIN. Set `tropo_method='a_linear'`. Save, Start.
5. Expected:
   - No degradation warning.
   - `aps_linear` runs (console output contains "loading the data").
   - Output contains `Atmosphere_a_linear_AOI_PS.mat` and `Atmosphere_a_linear_*.csv`.
   - Velocity values differ from a run with TRAIN unchecked.

### 3. GACOS correction (`a_gacos`) — optional, requires gacos.net data request

1. Same TRAIN install as above.
2. Launch `PHASE_StaMPS.mlapp`. Tick TRAIN. Set `tropo_method='a_gacos'`. Save, Start.
3. Expected:
   - MATLAB pauses at `keyboard` after creating `GACOS/` and `GACOS_download_info.txt`.
   - Download `.tar.gz` files from gacos.net per the printed instructions, place in `GACOS/`, press Continue in MATLAB editor.
   - MATLAB extracts/distributes `.ztd` files.
   - Output contains `Atmosphere_a_gacos_AOI_PS.mat` and `Atmosphere_a_gacos_*.csv`.

## Updates
- *April 2026*: Added interactive geographic map GUI for automatic AOI sub-setting. Introduced meteorologically-aware master image selection using Open-Meteo API. Automated parameter metadata detection for StaMPS. Dropped legacy Python 2.7 support.
- *March 2026*: Introduced Module 2 for geospatial PSI data analysis with deterministic and stochastic modeling.
- *September 2024*: Added *macOS* compatibility to the preprocessing application and improved master error handling.
<img width="2100" height="1181" alt="GitHubUpdates" src="https://github.com/user-attachments/assets/376d8baf-72b6-42d1-a08c-9d50f764bb28" />

## Planned updates
- Improve border constraints based on user selection.
- Introduce the handling of jumps in the displacement models.
- Introduce the possibility for a NRT processing.
- Native integration of Alaska Vertex API for direct Sentinel-1 image downloads.

## Acknowledgments
Special thanks to Jose Manuel Delgado Blasco and Dr. Michael Foumelis for the snap2stamps[^1] tool, and Prof. Andy Hooper for the StaMPS[^2] development. <br>

When using this software, please refer to:<br>
Monti, R., & Rossi, L. (2025). PHASE: a Matlab-based software for the DInSAR PS processing. Geodesy and Cartography, 51(2), 88–99. https://doi.org/10.3846/gac.2025.21995

[^1]: Foumelis, M., Delgado Blasco, J. M., Desnos, Y. L., Engdahl, M., Fernández, D., Veci, L. Lu, J. and Wong,
C. “SNAP - StaMPS Integrated processing for Sentinel-1 Persistent Scatterer Interferometry”. In
Geoscience and Remote Sensing Symposium (IGARSS), 2018 IEEE International, IEEE. <br>
[^2]: Hooper, A., A multi-temporal InSAR method incorporating both persistent scatterer and small baseline approaches, Geophys. Res. Lett., 35, L16,302, doi:10.1029/2008GL03465, 2008.
