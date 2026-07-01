"""One-shot script surfacing the silent TRAIN-degradation reset in the GUI.

Context
-------
When the user asks for a tropospheric correction (subtr_tropo='y') but TRAIN
is not on the MATLAB path, the Windows-port guard (Change #1, see
apply_train_windows_port.py) degrades `train_flag` to 1. Execution then takes
the `else` branch of `if train_flag == 0`, whose first statement is
`setparm('subtr_tropo', 'n');` — this is the second, overriding definition of
subtr_tropo the user sees in the Command Window ("y" then "n").

That override is intentional (no TRAIN => nothing to subtract), but until now
it was ONLY announced via `warning()` in the Command Window; the app's message
area said nothing, so users read the y->n flip as a bug (support ticket from
G. Elli, thesis under R. Monti).

This edit adds a single `updateOutput(app, ...)` notice at the reset site,
gated on `~train_available && subtr_tropo=='y'`, so the GUI explains WHY the
correction was skipped and HOW to enable it. Behavior is unchanged; the notice
is purely informational.

Deliberately NOT placed inside the `% begin/end TRAIN availability check`
sentinels: tests/test_phase_train_runtime_degradation.py extracts that block
and runs it standalone (no `app` object), so an `updateOutput(app, ...)` there
would break it. The notice lives ~180 lines later, at the consumption site.

Run from anywhere: python tools/apply_train_visibility_fix.py
Idempotent: re-runs detect the already-applied edit and exit 0.
"""
import sys
from pathlib import Path

TOOLS_DIR = Path(__file__).parent
REPO_ROOT = TOOLS_DIR.parent
sys.path.insert(0, str(TOOLS_DIR))
from mlapp_roundtrip import edit_mlapp

MLAPP = REPO_ROOT / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"
if not MLAPP.exists():
    raise SystemExit(
        f"Expected mlapp at {MLAPP}; is the script in the correct tools/ dir? "
        f"TOOLS_DIR={TOOLS_DIR}, REPO_ROOT={REPO_ROOT}"
    )

# Anchor: the `else` (train_flag != 0) branch and its first two statements.
# `setparm('subtr_tropo', 'n');` occurs exactly once in the file, which makes
# this block unique. Indentation is ground-truth: else=20 spaces, body=24.
BEFORE = (
    "                    else\n"
    "                        setparm('subtr_tropo', 'n');\n"
    "                        stamps(stamps_first_step,7); % execute StaMPS steps from first step to 7\n"
)

AFTER = (
    "                    else\n"
    "                        % begin TRAIN-degradation GUI notice (surfaces the otherwise-silent\n"
    "                        % subtr_tropo='n' reset; mirrors the Command Window warning above so\n"
    "                        % the user is not left thinking the requested correction actually ran).\n"
    "                        if ~train_available && strcmpi(subtr_tropo, 'y')\n"
    "                            updateOutput(app, ['NOTE: TRAIN was not found on the MATLAB path, so the requested ', ...\n"
    "                                'tropospheric correction (subtr_tropo=y, tropo_method=', tropo_method, ') is being ', ...\n"
    "                                'SKIPPED and subtr_tropo is forced to ''n'' for this run. To enable it, run ', ...\n"
    "                                'which(''aps_weather_model'') in MATLAB: if it is empty, add TRAIN with ', ...\n"
    "                                'addpath(genpath(''<TRAIN>/matlab'')); savepath, then press Start again. ', ...\n"
    "                                'For tropo_method=a_gacos also install GMT and make sure ''gmt --version'' works.']);\n"
    "                        end\n"
    "                        % end TRAIN-degradation GUI notice\n"
    "                        setparm('subtr_tropo', 'n');\n"
    "                        stamps(stamps_first_step,7); % execute StaMPS steps from first step to 7\n"
)


def edit(xml: str) -> str:
    if "% begin TRAIN-degradation GUI notice" in xml:
        print("  visibility notice: already applied, skipping")
        return xml
    n = xml.count(BEFORE)
    if n != 1:
        raise SystemExit(
            f"  visibility notice: BEFORE block found {n} times; expected 1. "
            f"Aborting (file structure may have changed)."
        )
    print("  visibility notice: applying")
    return xml.replace(BEFORE, AFTER, 1)


if __name__ == "__main__":
    edit_mlapp(MLAPP, edit)
    print("Done. Verify with: python tools/apply_train_visibility_fix.py  (should say 'already applied')")
