"""Repair the PHASE_StaMPS App Designer regression introduced in 290e56c.

Commit 290e56c was saved from a stale App Designer working copy.  Besides its
intended TRAIN/subtr_tropo compatibility check, that save removed the Windows
StaMPS branches and several other changes already present on main.  Restore the
last complete app from the merge immediately before the regression, then apply
the intended compatibility check as a narrow source edit.

The script is deliberately guarded: it only replaces an app that has the known
regression signature, and it validates the historical archive before writing.
Run from anywhere with:

    python tools/repair_phase_stamps_regression.py
"""

from __future__ import annotations

import io
import subprocess
import sys
import zipfile
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TOOLS_DIR.parent
sys.path.insert(0, str(TOOLS_DIR))

from mlapp_roundtrip import edit_mlapp


MLAPP_REL = Path("PHASE_Preprocessing/PHASE_StaMPS.mlapp")
MLAPP = REPO_ROOT / MLAPP_REL
GOOD_REVISION = "918a4caa481ebdd747c3e44caba7f574ff4d2059"

REGRESSION_MARKERS = (
    "AutoDetectParameters(app);is",
    "stamps_path = which('StaMPS_CONFIG.bash');",
    "system(strjoin({source_stamps, source_train, source_snap}, ';'));",
)

RESTORED_MARKERS = (
    "StaMPS_CONFIG.ps1",
    "mt_prep_snap.bat",
    "% begin TRAIN availability check",
    "% begin TRAIN interactive recovery",
    "% begin TRAIN-degradation GUI notice",
    "self-bootstrapping .bat shim",
)

SAFETY_BEFORE = (
    "                if contains(ph_output, 'unwrapped')\n"
    "                    if train_flag == 0\n"
)
SAFETY_AFTER = (
    "                if contains(ph_output, 'unwrapped')\n"
    "                    % TRAIN is usable only when both controls request correction.\n"
    "                    if train_flag == 0 && strcmpi(strtrim(subtr_tropo), 'y')\n"
)


def document_xml(archive: bytes) -> str:
    with zipfile.ZipFile(io.BytesIO(archive)) as mlapp_zip:
        return mlapp_zip.read("matlab/document.xml").decode("utf-8")


def apply_safety_check(xml: str) -> str:
    if SAFETY_AFTER in xml:
        return xml
    count = xml.count(SAFETY_BEFORE)
    if count != 1:
        raise SystemExit(
            f"Safety-check anchor found {count} times; expected exactly 1. Aborting."
        )
    return xml.replace(SAFETY_BEFORE, SAFETY_AFTER, 1)


def main() -> None:
    current = MLAPP.read_bytes()
    current_xml = document_xml(current)

    if all(marker in current_xml for marker in RESTORED_MARKERS):
        if SAFETY_AFTER not in current_xml:
            edit_mlapp(MLAPP, apply_safety_check)
            print("App was already restored; applied the compatibility check.")
        else:
            print("PHASE_StaMPS.mlapp is already repaired; nothing to do.")
        return

    missing_regression_markers = [
        marker for marker in REGRESSION_MARKERS if marker not in current_xml
    ]
    if missing_regression_markers:
        raise SystemExit(
            "The app does not match the known regression signature; refusing to "
            "replace it. Missing markers: " + repr(missing_regression_markers)
        )

    result = subprocess.run(
        ["git", "show", f"{GOOD_REVISION}:{MLAPP_REL.as_posix()}"],
        cwd=REPO_ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    restored_xml = document_xml(result.stdout)
    missing_restored_markers = [
        marker for marker in RESTORED_MARKERS if marker not in restored_xml
    ]
    if missing_restored_markers:
        raise SystemExit(
            "Historical app failed validation; refusing to replace the current app. "
            "Missing markers: " + repr(missing_restored_markers)
        )

    temporary = MLAPP.with_suffix(".mlapp.restore-tmp")
    temporary.write_bytes(result.stdout)
    temporary.replace(MLAPP)
    edit_mlapp(MLAPP, apply_safety_check)
    print(
        "Restored the complete pre-regression app and applied the "
        "TRAIN/subtr_tropo compatibility check."
    )


if __name__ == "__main__":
    main()
