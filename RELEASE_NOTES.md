# PHASE v6.0

PHASE v6.0 promotes the new standalone applications for preprocessing,
StaMPS processing and geospatial modelling to the main production release.
The complete workflow has been validated end-to-end on Windows with real data.

## Windows installation — important

**Windows users should download and run `install-phase.exe` attached to this
release. Do not clone the repository for a normal Windows installation.**

The installer creates the complete working environment and the three PHASE
shortcuts. It installs or configures PHASE, StaMPS, TRAIN, the required native
Windows executables, Python and the MATLAB paths. A plain repository clone does
not perform these operations and is intended only for PHASE development.

After installation, launch the applications using:

- `PHASE Preprocessing.lnk`
- `PHASE StaMPS.lnk`
- `PHASE Model.lnk`

The editable MATLAB source remains available in the visible `engine` folder.

## Highlights

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
