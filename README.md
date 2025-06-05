# PHASE (Persistent scatterer Highly Automated Suite for Environmental monitoring)

**PHASE** (**P**ersistent scatterer **H**ighly **A**utomated **S**uite for **E**nvironmental monitoring) is a MATLAB-based software suite for the automatization of the InSAR PSI processing. <br>
It is based on the well-known and widely used *snap2stamps* and *StaMPS* software, while coming with a set of new features and improvements.

![Logo](https://github.com/user-attachments/assets/1ced30e1-e8c5-4186-bc7f-4c6dec5ca12d)

## SAR Satellites compatibility
- Sentinel-1 (from European Space Agency)
- COSMO-SkyMed (from Agenzia Spaziale Italiana)

## Required Software
- SNAP
- MATLAB
- StaMPS

## Required OS
- *Linux*: it is mandatory as StaMPS only works in this environment. The entire processing can be executed from start to end in it.
- *Windows*: only the preprocessing part can be executed in it.
- *macOS*: only the preprocessing part can be executed in it.

## Installation and Setup

> [!NOTE]
> A more detailed step-by-step guide is available in the provided user manual. <br>
> Before start using PHASE please read carefully the whole manual!

### Preliminary Steps
1. **Install SNAP Software** <br>
   Download and install [SNAP](https://step.esa.int/main/download/snap-download/) from the European Space Agency website. <br> <br>
   Check that all the mandatory SNAP plugin modules are installed:
   - Microwave Toolbox Kit Module
   - Optical Toolbox Kit Module
   - SMOS-Box Kit Module
   - Radarsat Polarimetric Toolkit Module
   - ESA SNAPPY
   - EOMTBX<br>

    After installing SNAP on your computer, you are suggested to review the parameters set in:
     - $HOME/snap/bin/gpt.vmoptions and modify the param
          - –Xmx 12G (according to your computer set up; i.e –Xmx 512M )
     - $HOME/snap/etc/snap.properties
          - #snap.home=
          - #snap.userdir=
          - snap.jai.tileCacheSize = 1024
          - snap.jai.defaultTileSize = 512

3. **Install Required Python Modules:** <br>
   Install [Python](https://www.python.org/downloads/) 2.7 or 3.x on your machine. Ensure you have the necessary libraries; if not, install them with:
   ```
   pip install os
   pip install sys
   pip install pathlib
   pip install shutil
   pip install glob
   pip install subprocess
   pip install shlex
   pip install time
   ```
4. **Install xterm (only Linux Users):** <br>
   Install xterm by running `sudo apt-get install xterm` in the terminal.

5. **Install StaMPS:** <br>
   Install [StaMPS](https://homepages.see.leeds.ac.uk/~earahoo/stamps/) from the official GitHub repository.
   ```
   git clone https://github.com/dbekaert/StaMPS/releases/tag/v4.1-beta
   ```

6. **PHASE suite**
   - Download the *PHASE_MANUAL* and the *PHASE_python2* or *PHASE_python3* folder based on your python version.
   - Move or copy the downloaded folder (*PHASE_pythonX*) in your project folder, anywhere on your computer.
   - Exectute the *PS_InSAR_Preprocessing.mlapp* MATLAB application.
   - Tune all the configurable parameters across all the available tabs.
   - Once the preprocessing is complete, execute the *StaMPS_Automate.mlapp* MATLAB application (in *Linux* it will automatically open upon completion).

### Processing Steps

1.	**Automated SAR Images Download:** <br>
Download images from the Alaska SAR satellite facility for Sentinel-1 or manually from the Italian Space Agency website for COSMO-SkyMed.
2.	**Master Image Processing:** <br>
Split and apply orbit correction. All the configurable parameters can be set through the GUI.
3.	**Slaves Pre-Processing:** <br>
Steps include slaves preparation, splitting and orbit correction, coregistration and interferogram formation, StaMPS export, average scene intensity computation, and local incidence angle and coherence images computation. All the configurable parameters can be set through the GUI.
4.	**StaMPS Processing:** <br>
Data preparation, parameter definition, and execution of StaMPS PS analysis. Displacement time series export in Excel format.

## Possible Errors and Solutions
The procedure has been tested on SNAP 10.x, Python 2.7, Python 3.11, Ubuntu 20.04, Windows 10, macOS Sequoia (15.1), and MATLAB 2023b. <br>
> [!TIP]
> Refer to the manual for solutions to common errors encountered during the StaMPS processing.

## Updates
- *September 2024*: added *macOS* compatibility to the preprocessing application.
- *September 2024*: improved master error handling.
![GitHubUpdates](https://github.com/user-attachments/assets/7324e1e9-a191-4dff-bff8-03a27e8040e9)

## Planned updates
- add the support to SAOCOM data processing.
- next major release will introduce the spatio-temporal time series processing with deterministic and stochastic methods. The same OS compatibility will be kept also for this new module.

## Acknowledgments
Special thanks to Jose Manuel Delgado Blasco and Dr. Michael Foumelis for the snap2stamps[^1] tool, and Prof. Andy Hooper for the StaMPS[^2] development. <br>

When using this software, please refer to:<br>
Monti, R., & Rossi, L. (2025). PHASE: a Matlab-based software for the DInSAR PS processing. Geodesy and Cartography, 51(2), 88–99. https://doi.org/10.3846/gac.2025.21995

[^1]: Foumelis, M., Delgado Blasco, J. M., Desnos, Y. L., Engdahl, M., Fernández, D., Veci, L. Lu, J. and Wong,
C. “SNAP - StaMPS Integrated processing for Sentinel-1 Persistent Scatterer Interferometry”. In
Geoscience and Remote Sensing Symposium (IGARSS), 2018 IEEE International, IEEE. <br>
[^2]: Hooper, A., A multi-temporal InSAR method incorporating both persistent scatterer and small baseline approaches, Geophys. Res. Lett., 35, L16,302, doi:10.1029/2008GL03465, 2008.
