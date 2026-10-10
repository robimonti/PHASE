# PHASE 7 public-release checklist

The `v7.x.y` tag workflow builds a Windows EXE tied to that tag, a Linux
x86_64 AppImage installer, and the `phase7-engine.zip` in-app update asset.
It creates a **draft** GitHub release only. Nothing is published automatically.

## Before tagging

1. Merge the intended release commit to `main`; run the Python suite and
   inspect `git diff --check` and the Windows installer syntax check.
2. Confirm the Linux native-runtime smoke workflow passes on Ubuntu. This
   compiles StaMPS/Triangle and checks TRAIN, but **does not** run MATLAB,
   SNAP or a real Linux PSI stack. Label Linux accordingly in the release.
3. Review third-party source and binary licenses, especially Triangle and
   SNAPHU, before distributing any bundled runtime. MATLAB/SNAP remain
   external prerequisites and are not included in PHASE installers.
4. Build the full Apple Silicon DMG from the pinned runtime manifest. Include
   matching GNU awk and SNAPHU source archives. This release is intentionally
   unsigned/not notarized; document Gatekeeper's **Open Anyway** path and do
   not describe the DMG as frictionless or Apple-verified.
5. Ensure the version displayed in the hub, Windows EXE, macOS app/installer,
   Git tag and release notes agrees. Test a clean install and one update from
   a previous managed PHASE 7 installation.
6. Review the downloadable user manual: the current PDF predates the PHASE 7
   unified workflow and should be updated or clearly labeled historical.

## Draft-release review

1. Confirm `install-phase.exe`, `PHASE-7-linux-x86_64.AppImage` and
   `phase7-engine.zip` are attached to the draft release. Attach the complete
   `PHASE-7-macos-arm64.dmg` built locally with the complete runtime and source
   archives; the workflow does not manufacture a Mac runtime or publish an
   engine-only DMG.
2. Attach `SHA256SUMS.txt`, the matching GNU awk and SNAPHU source archives,
   verify SHA-256 of every uploaded asset, and confirm GitHub reports a
   `sha256:` digest for `phase7-engine.zip` (required by the in-app updater).
3. Download the release assets as a user would. Check that Windows installs
   the tagged PHASE engine, macOS launches from `~/Applications/PHASE.app`,
   Linux shows its prerequisites and builds the pinned StaMPS/TRAIN revisions,
   and each installation reports the release tag rather than `dev`.
4. Publish the draft only after the Mac asset, license review, version checks
   and update test are complete. If Linux PSI has not been tested, retain its
   explicit unvalidated/experimental status instead of claiming parity.

Existing PHASE 7 projects are **not** renamed. New projects use the sequential
`01_INPUT`, `02_PROCESSING_INTERNAL`, `03_RESULTS`, `04_LOGS` layout; the old
`IMPORTED_RESULTS` catalog lives under `04_LOGS` only for a legacy import.
