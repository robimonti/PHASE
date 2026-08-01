# PHASE Model 6.0 — implementation notes

`PHASE_Model.m` is the production entry point for Module 2. It wraps the
validated `PHASE_Model_beta.m` implementation, whose internal name is retained
for backward compatibility, and does not load `PHASE_model.mlapp` at runtime.

The validated App Designer class is mechanically extracted into
`+phase_model_beta/LegacyEngine.m`, preserving all callbacks and the complete
scientific workflow while making the source inspectable and versionable.
`+phase_model_beta/App.m` keeps that engine hidden and presents the same
white, minimal Liquid Glass shell used by the preprocessing and StaMPS modules.
The visible configuration, calendar, conditional controls, save/load state,
progress and live diary are implemented in `phase_model_beta_ui`.

The AOI page includes the same satellite map used by preprocessing. It can
estimate the complete PS extent or accept an arbitrary polygon drawn and
refined directly on the map; the four visible coordinates are its bounding
box, while the complete polygon is retained for PS selection. Model-specific pages are shown only when
they apply to the processing family and project dimension. Temporal workflows
also expose Automatic/Manual controls for the five scientific thresholds in
`MatlabFunctions/ModellingInTime.m`; Automatic keeps the function defaults,
while Manual passes the visible value to the processing backend.

Regenerate it after an intentional stable-model change with:

```text
python tools/extract_model_beta_engine.py
```

Launch from the PHASE root:

```matlab
addpath(genpath(pwd))
phase_model_beta.selfTest
PHASE_Model
```

On Windows, the first launch automatically reads the Python interpreter
selected by the PHASE installer from `%APPDATA%\PHASE\python.txt`. A Python
path saved in `input_model.mat` remains authoritative.

The legacy waitbar and report-formatting terminal windows are not used by the
production application. Progress and MATLAB/Python messages remain inside **Run monitor**. The
monitor also includes **Hard stop now**, which interrupts external helpers and
stops MATLAB at the next safe numerical checkpoint. Pure temporal processing
does not allocate a spatial grid, and each report figure is verified on disk
before it is embedded in Excel. App Designer callbacks are preserved in the
extracted engine; the editable
scientific functions remain normal MATLAB files under `MatlabFunctions`.
