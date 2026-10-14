# PHASE Model 6.0 — implementation notes

`PHASE_Model.m` is the production entry point for Module 2. It wraps the
validated `PHASE_Model_beta.m` implementation, whose internal name is retained
for backward compatibility, and does not load `PHASE_model.mlapp` at runtime.

The processing runner is now the editable headless MATLAB class
`+phase_model_beta/ProcessingEngine.m`. It does not instantiate an App Designer
window or read an MLAPP at runtime. `+phase_model_beta/App.m` presents the same
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

Change the `.m` runner and its `MatlabFunctions` dependencies directly. The
former App Designer source is archived only for historical reference.

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
before it is embedded in Excel. The scientific functions remain normal MATLAB
files under `MatlabFunctions`.
