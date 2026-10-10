# PHASE

**Persistent Scatterer Highly Automated Suite for Environmental Monitoring**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19605360.svg)](https://doi.org/10.5281/zenodo.19605360)
[![Latest release](https://img.shields.io/github/v/release/robimonti/PHASE)](https://github.com/robimonti/PHASE/releases/latest)

PHASE is a MATLAB-based application for InSAR persistent-scatterer processing
and displacement analysis. Its unified hub brings **Preprocessing**, **StaMPS
PSI** and **Displacement Modeling** into one window.

> [!IMPORTANT]
> Install from [GitHub Releases](https://github.com/robimonti/PHASE/releases).
> Do not use **Code → Download ZIP** as an installer. The macOS installer is not Apple Developer-ID
> signed or notarized; see the opening instructions below. Linux remains
> experimental until a complete PSI run is tested there.

## Install

Download the installer matching your OS and version from
[GitHub Releases](https://github.com/robimonti/PHASE/releases):

| OS | PHASE 7 installer | Validation status |
| --- | --- | --- |
| Windows 10/11 | `install-phase.exe` | Guided workflow tested |
| macOS Apple Silicon | `PHASE-7-macos-arm64.dmg` | Preprocessing and StaMPS PSI tested, including Step 8 and TRAIN `a_linear` |
| Linux x86_64 | `PHASE-7-linux-x86_64.AppImage` | Native build path prepared; end-to-end PSI not yet tested |

Apple Intel is outside the supported scope. MATLAB and ESA SNAP are external
prerequisites. The Linux installer also needs a C++ build toolchain, CMake,
`snaphu`, `gawk`, `csh`, and `zenity` or `kdialog`; it builds the pinned StaMPS
and TRAIN forks locally. See [installation requirements](installer/README.md).

On macOS, open the downloaded DMG and copy/open **PHASE Installer.app**. If
macOS blocks this unsigned app, try to open it once, then use **System Settings
→ Privacy & Security → Open Anyway** and confirm **Open**. Only do this for
an asset downloaded from the PHASE GitHub release; verify its published
SHA-256 first. Do not disable Gatekeeper globally. The unsigned DMG is a
trade-off made to avoid a paid Apple Developer membership, not a notarized
package. The third-party source archives and license notices are supplied
with the DMG and release downloads.

Install PHASE once per computer, then create as many projects as needed in
folders of your choice. The hub is available as **PHASE** on Windows,
`~/Applications/PHASE.app` on macOS, and in the application menu on Linux.

## Project folders

New projects have four top-level folders:

```text
Project/
  01_INPUT/                 Optional SAR archive and AOI files
  02_PROCESSING_INTERNAL/   PHASE working files — do not edit during a run
  03_RESULTS/               Final PSI series, modeling, figures, GIS, reports
  04_LOGS/                  Run logs and optional legacy-import catalog
```

Project files stay separate from the installed PHASE engine. Older PHASE 7
projects keep their original folder names; PHASE does not rename paths that
existing processing state may reference. The former `IMPORTED_RESULTS` folder
only held a catalog of results imported from legacy workspaces; new projects
store that catalog in `04_LOGS` only when an import is performed. Details are
in the [project architecture](Documenti/PHASE_7_Architettura_Progetti.md).

## Updates

**Check for updates** in the hub downloads and verifies a compatible stable
v7 engine from GitHub. Close PHASE and launch it again from the installed
shortcut to apply the update. Project folders and installed StaMPS/TRAIN
runtimes are left untouched. The first compatible update asset will arrive
with the first v7 release. The previous engine is retained as a backup.

## Science and documentation

PHASE supports Sentinel-1 and COSMO-SkyMed processing. The PSI workflow uses
[ESA SNAP](https://step.esa.int/main/download/snap-download/) and
[StaMPS](https://github.com/pyccino/StaMPS), with optional
[TRAIN](https://github.com/pyccino/TRAIN) atmospheric correction. An
end-to-end software run establishes operational compatibility, not scientific
accuracy for every dataset; users should inspect warnings and validate
deformation estimates against independent observations where possible.

- [PHASE 6 manual (historical; PHASE 7 revision pending)](PHASE_Manual.pdf)
- [Release notes](RELEASE_NOTES.md)
- [Installer and release preparation](installer/README.md)
- [Scientific export corrections](docs/PSI_EXPORT_CORRECTIONS.md)
- [Third-party components and source code](Documenti/PHASE_7_THIRD_PARTY.md)
- [Archived detailed PHASE 6 guide](Documenti/PHASE6_GUIDE.md)

The archived standalone PHASE 6 applications are in [`legacy/`](legacy/README.md)
for provenance; new users should open the unified PHASE hub.
