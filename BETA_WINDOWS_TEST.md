# PHASE standalone beta — Windows acceptance test

This package contains the three standalone MATLAB applications and their
editable backends. It intentionally excludes the legacy `.mlapp` files.

## 1. Clean test folder

Extract the archive into a new local folder, for example:

```text
E:\PHASE_beta_complete_test
```

Do not copy only the launcher `.m` files: the package folders beside them are
part of the runtime.

## 2. Prepare StaMPS and TRAIN

The standalone transfer archive does not embed third-party git repositories.
Before real StaMPS processing, open PowerShell in the extracted package root
and run:

```powershell
powershell -ExecutionPolicy Bypass -File .\prepare-beta-test-windows.ps1
```

The script reproduces the relevant installer work without touching datasets:
it clones/updates `StaMPS` and `TRAIN` beside PHASE and installs the nine
mandatory StaMPS Windows executables. PHASE StaMPS automatically detects this
layout and adds it to the current MATLAB process; a permanent `savepath` is
not required for the acceptance test.

## 3. Smoke tests

Open MATLAB, then run:

```matlab
cd('E:\PHASE_beta_complete_test')
addpath(genpath(pwd))

phase_preprocessing_beta.selfTest
phase_stamps_beta.selfTest
phase_model_beta.selfTest
```

All three results must report `ok: 1`.

## 4. Preprocessing and StaMPS

Run:

```matlab
PHASE_Preprocessing_beta
```

Use a disposable copy of a known dataset. Verify that preprocessing output is
streamed into **Run monitor** and that no separate CMD window opens. At the
StaMPS handoff, verify:

- the selected first/last step range is drawn from both ends in the sidebar;
- a manually selected master date uses the calendar picker;
- StaMPS/third-party MATLAB messages appear incrementally in **Run monitor**;
- PHASE-owned commands do not open CMD windows;
- **TS Points** opens inside the same PHASE window and **Back to PHASE**
  returns to the main interface.
- Step 3 displays **Random percentage** for `PERCENT` and **Random density**
  for `DENSITY`; StaMPS receives the corresponding parameter.

If StaMPS or one of its native executables is missing, Start must now stop
before scientific processing and display the exact missing dependency in
PHASE's Run monitor/error dialog.

## 5. PHASE Model

From the package root run:

```matlab
PHASE_Model_beta
```

The complete processing implementation is readable in
`+phase_model_beta\LegacyEngine.m`; the launcher does not load a legacy app.
The backend stays hidden behind the same modern shell as Modules 1A and 1B.
Verify that its calendar, conditional parameter cards and **Run monitor**
remain inside the PHASE Model window and that report formatting opens no CMD
window.
Also verify:

- **Save** creates/updates `input_model.mat` without an input-argument error;
- the bounding-box selector remains full width and its four coordinates form
  a 2 × 2 longitude/latitude layout;
- **Select full PS extent** estimates all four bounds from the selected input;
- the AOI map loads its satellite background, accepts a non-rectangular
  polygon, and updates the four bounding coordinates without discarding the
  polygon shape;
- changing **Processing family** and **Project dimension** immediately shows
  only the applicable Observations, Temporal, 1D or 2D parameter pages;
- every Automatic/Manual temporal threshold shows its value only in Manual
  mode;
- the first Output option reads **Export interpolated PS observations**.
- **Hard stop now** is enabled only during a run and stops at the next safe
  numerical checkpoint;
- a temporal-only run creates and embeds the General, Raw displacement,
  Modelled displacement and Uncertainty figures in the Excel report without
  allocating the 1D/2D interpolation grid.

On an installer-created system it reads the selected Python path from
`%APPDATA%\PHASE\python.txt`, unless `input_model.mat` already contains an
explicit value.

## 6. Installer source

`installer\install-phase.ps1` is the beta installer source. It clones
`codex/phase-stamps-beta`, validates that the complete standalone runtime is
present, removes legacy/development files from the installed engine, creates
the three visible MATLAB launch shortcuts and leaves the editable engine folder visible. The
StaMPS shortcut asks for an explicit `ASC_*`/`DSC_*` dataset folder.

Do not use the old `installer\install-phase.exe` from another checkout: an
updated EXE must be compiled on Windows after the tested beta branch has been
pushed:

```powershell
powershell -ExecutionPolicy Bypass -File installer\compile-to-exe.ps1
```
