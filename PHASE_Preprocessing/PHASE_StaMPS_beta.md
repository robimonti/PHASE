# PHASE StaMPS 6.0 — implementation notes

`PHASE_StaMPS` is the production entry point. It wraps the validated
`PHASE_StaMPS_beta` implementation, whose internal name is retained for
backward compatibility. It uses the same `input_StaMPS.mat`, StaMPS folders
and output products as the former app, but does not load or copy it at runtime.

## Launch

In MATLAB:

```matlab
addpath('C:\path\to\PHASE\PHASE_Preprocessing')
PHASE_StaMPS('C:\path\to\ASC_Jul25_May26')
```

Before the first launch, the configuration layer can be checked without
running StaMPS:

```matlab
phase_stamps_beta.selfTest
```

If MATLAB's current folder already is the `ASC_*`/`DSC_*` processing folder:

```matlab
PHASE_StaMPS
```

When no `input_StaMPS.mat` exists in the selected folder, the interface opens
with an initial configuration. It detects the PHASE project and acquisition
metadata when available; select the StaMPS installation folder, review the
values and press **Save** to create `input_StaMPS.mat`.

Copying only `PHASE_StaMPS_beta.m` into an `ASC_*`/`DSC_*` folder is optional:
the copied launcher resolves the canonical editable runtime from the sibling
`PHASE_Preprocessing` folder. The full `+phase_stamps_beta` and
`phase_stamps_beta_ui` folders must remain in the PHASE installation.

## Safety model

- The former `PHASE_StaMPS.mlapp` is archived under `legacy` for provenance.
- The standalone engine reads and writes the same MAT variables.
- Saving is transactional: the previous MAT is copied to
  `input_StaMPS.mat.bak` before replacement.
- **Start processing** is disabled while visible values are unsaved.
- MATLAB validates paths, dates, temporal windows and the selected StaMPS
  range again before execution.
- On Windows, PHASE-owned external commands run without visible CMD windows.
  MATLAB/StaMPS output is mirrored incrementally into **Run monitor** through
  a per-session diary, including messages printed below the PHASE layer.
- The master date uses a native calendar control. The selected first/last
  StaMPS range is represented from both ends in the sidebar progress bar.
- The first standalone engine was mechanically extracted from the validated
  Start callback. A regression test proves that all `setparm`, `setparm_aps`
  and `stamps` calls remain identical.

## Source layout

```text
PHASE_Preprocessing/
├── PHASE_StaMPS_beta.m                 launcher
├── +phase_stamps_beta/
│   ├── App.m                           MATLAB/HTML controller
│   ├── defaultConfig.m                 processing defaults
│   ├── schema.m                        labels, groups, units and options
│   ├── loadConfig.m / saveConfig.m     MAT compatibility
│   ├── configToUi.m / uiToConfig.m     typed UI conversion
│   ├── validateConfig.m                safety checks
│   ├── autoDetectConfig.m              SNAP metadata detection
│   ├── LegacyAppAdapter.m              stable-engine compatibility
│   ├── runProcessing.m                 editable StaMPS engine
│   ├── runCommandHidden.m               no-console streamed commands
│   └── openTsPicker.m                  TS Points integration
└── phase_stamps_beta_ui/
    ├── index.html
    ├── styles.css
    └── app.js
```

There is no production `.mlapp` and no `document.xml` patch. Labels and parameter
metadata can be changed in `schema.m`; visual styling lives in `styles.css`;
interactions live in `app.js`; processing remains in MATLAB. The archived
`PHASE_StaMPS.mlapp` is not loaded, copied or otherwise required at runtime.

## Windows acceptance test

Use a copy of an already completed `ASC_*`/`DSC_*` folder.

1. Launch PHASE StaMPS with the explicit folder path.
2. Confirm that all saved values load correctly, including the three temporal
   windows and first/last StaMPS steps.
3. Change one harmless value, confirm the amber **Unsaved** state, then press
   **Load** and verify that the saved value returns.
4. Change it again, press **Save**, close/reopen the app and verify persistence.
5. Without TRAIN, run Steps 6 → 8 and confirm the sidebar bar begins at Step 6,
   the live log reports the selected range and completes Step 8 without
   opening CMD windows.
6. Repeat the production workflow with TRAIN/GACOS enabled.
7. Open **TS Points** and confirm the picker replaces the central workspace
   inside the same PHASE window; **Back to PHASE** must restore the main UI.
8. Compare the generated CSV/XLSX and principal MAT products with the stable
   app on the same input dataset.

The complete release acceptance procedure is documented in
`WINDOWS_ACCEPTANCE_TEST.md`.
