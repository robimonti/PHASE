"""One-shot script applying the three edits for the TRAIN Windows port.

Run from anywhere: python tools/apply_train_windows_port.py
Path-anchored on __file__ so cwd doesn't matter.
Idempotent: re-runs detect already-applied edits and exit 0.
"""
import sys
from pathlib import Path
TOOLS_DIR = Path(__file__).parent
REPO_ROOT = TOOLS_DIR.parent
sys.path.insert(0, str(TOOLS_DIR))
from mlapp_roundtrip import edit_mlapp

MLAPP = REPO_ROOT / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"
# Use an explicit SystemExit rather than `assert` so `python -O` cannot
# strip the check.
if not MLAPP.exists():
    raise SystemExit(
        f"Expected mlapp at {MLAPP}; is the script in the correct tools/ dir? "
        f"TOOLS_DIR={TOOLS_DIR}, REPO_ROOT={REPO_ROOT}"
    )

# Change #1 — insert guard after input_StaMPS.mat load()
CHANGE1_BEFORE = (
    "                    'ref_radius', 'ref_velocity', 'plot_s', "
    "'ref_centre_lonlat_w', 'ref_radius_w', 'ph_output');\n"
    "\n"
    "                % Conversion of variables' format\n"
)
CHANGE1_AFTER = (
    "                    'ref_radius', 'ref_velocity', 'plot_s', "
    "'ref_centre_lonlat_w', 'ref_radius_w', 'ph_output');\n"
    "\n"
    "                % begin TRAIN availability check (Windows port, cross-platform no-op\n"
    "                % when TRAIN is correctly installed). Multi-function probe\n"
    "                % avoids false positives from unrelated aps_linear.m on path\n"
    "                % or incomplete installs. On miss, degrade silently using a\n"
    "                % LOCAL train_flag override — do not touch app.train_flag nor\n"
    "                % the checkbox, so the user's original intent is preserved\n"
    "                % for the next run (they can install TRAIN and press Start\n"
    "                % again without re-ticking the box).\n"
    "                train_available = ~isempty(which('aps_linear')) && ...\n"
    "                                  ~isempty(which('aps_weather_model')) && ...\n"
    "                                  ~isempty(which('setparm_aps'));\n"
    "                if train_flag == 0 && ~train_available\n"
    "                    warning('StaMPS:phase:trainNotAvailable', '%s', ...\n"
    "                        ['TRAIN not found on MATLAB path (or install is incomplete). ', ...\n"
    "                         'Install TRAIN (e.g. git clone https://github.com/dbekaert/TRAIN.git) ', ...\n"
    "                         'then in MATLAB: addpath(genpath(''<TRAIN>/matlab'')); savepath. ', ...\n"
    "                         'Proceeding without atmospheric correction.']);\n"
    "                    train_flag = 1;\n"
    "                end\n"
    "                % end TRAIN availability check\n"
    "\n"
    "                % Conversion of variables' format\n"
)

# Change #2 — remove misleading stub warning
CHANGE2_BEFORE = (
    "                        if isunix\n"
    "    system(strjoin({source_stamps, source_train}, ';'));\n"
    "else\n"
    "    % TRAIN Linux-only; on Windows just emit a warning\n"
    "    warning('StaMPS:phase:trainUnavailable', 'TRAIN is not available on Windows');\n"
    "end % source all the softwares and prepare the data\n"
)
CHANGE2_AFTER = (
    "                        if isunix\n"
    "    system(strjoin({source_stamps, source_train}, ';'));\n"
    "else\n"
    "    % Windows: when train_flag==0 reaches here, Change #1 has already\n"
    "    % verified TRAIN is on MATLABPATH; no shell config to source.\n"
    "end % source all the softwares and prepare the data\n"
)

# Change #3 — update stale comment
# Note: 4 leading spaces, NOT 20. The statement sits inside an unindented
# else-clause (`else\n    train_path = '';`) in the original file.
CHANGE3_BEFORE = (
    "    train_path = '';   "
    "% TRAIN not ported to Windows; PHASE skips if empty\n"
)
CHANGE3_AFTER = (
    "    train_path = '';   "
    "% Windows: TRAIN on MATLABPATH, no shell config to source\n"
)


def apply_edit(xml: str, before: str, after: str, name: str) -> str:
    if after in xml and before not in xml:
        print(f"  {name}: already applied, skipping")
        return xml
    if xml.count(before) != 1:
        raise SystemExit(
            f"  {name}: BEFORE block found {xml.count(before)} times; "
            f"expected exactly 1. Aborting."
        )
    print(f"  {name}: applying")
    return xml.replace(before, after, 1)


def edit(xml: str) -> str:
    xml = apply_edit(xml, CHANGE1_BEFORE, CHANGE1_AFTER, "Change #1")
    xml = apply_edit(xml, CHANGE2_BEFORE, CHANGE2_AFTER, "Change #2")
    xml = apply_edit(xml, CHANGE3_BEFORE, CHANGE3_AFTER, "Change #3")
    return xml


if __name__ == "__main__":
    edit_mlapp(MLAPP, edit)
    print("Done. Verify with: pytest tests/test_phase_train_detection.py")
