"""One-shot script adding an interactive "locate TRAIN" recovery dialog.

When the user asked for the tropospheric correction (train_flag == 0) but TRAIN
is not on the MATLAB path, the Windows-port guard degrades to train_flag = 1 and
the run proceeds WITHOUT correction (subtr_tropo forced to 'n'). That is safe but
opaque: the user just sees the correction vanish.

This edit adds, right AFTER the guard, an optional GUI recovery: a folder picker
(uigetdir) lets the user point at their TRAIN install; we addpath+validate it,
and on success re-enable the correction for this run (train_flag = 0) and offer
to make it permanent (savepath, with a userpath startup.m fallback when pathdef.m
is not writable). On cancel/skip the behaviour is unchanged (degrade).

Placement rules (do NOT break the tested guard):
- `train_requested` is captured BEFORE the `% begin TRAIN availability check`
  sentinel (the guard mutates train_flag, so we must snapshot the intent first).
- The interactive block sits AFTER the `% end TRAIN availability check` sentinel,
  so tests/test_phase_train_runtime_degradation.py (which extracts and runs the
  guard block standalone, with no `app`) is unaffected.

Run from anywhere: python tools/apply_train_locate_dialog.py
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
    raise SystemExit(f"Expected mlapp at {MLAPP}")

# --- Edit A: snapshot the user's intent just before the guard mutates train_flag.
BEGIN_SENTINEL = (
    "                % begin TRAIN availability check (Windows port, cross-platform no-op\n"
)
BEGIN_SENTINEL_NEW = (
    "                train_requested = (train_flag == 0);\n"
    + BEGIN_SENTINEL
)

# --- Edit B: the interactive recovery, inserted right after the end sentinel.
END_SENTINEL = "                % end TRAIN availability check\n"

# raw string: backslashes ('\n', '\') must reach document.xml verbatim so MATLAB
# fprintf/strrep see them; real newlines below are real newlines in the code.
RECOVERY = r'''
                % begin TRAIN interactive recovery (Windows-port UX). If the user
                % asked for the correction (train_flag was 0) but TRAIN is not on
                % the path, offer a folder picker to locate it and fix the run on
                % the fly, instead of silently proceeding without correction.
                if train_requested && ~train_available
                    reco_choice = uiconfirm(app.UIFigure, ...
                        ['TRAIN was not found on the MATLAB path, so the tropospheric correction (' tropo_method ') cannot run. ', ...
                         'Locate your TRAIN folder now to fix it, or skip the correction for this run.'], ...
                        'TRAIN not found', ...
                        'Options', {'Locate TRAIN folder...', 'Skip correction'}, ...
                        'DefaultOption', 1, 'CancelOption', 2);
                    if strcmp(reco_choice, 'Locate TRAIN folder...')
                        while true
                            train_dir = uigetdir('', 'Select your TRAIN folder (the one containing the "matlab" subfolder)');
                            if isequal(train_dir, 0)
                                updateOutput(app, 'TRAIN selection cancelled - proceeding without atmospheric correction.');
                                break
                            end
                            if exist(fullfile(train_dir, 'matlab'), 'dir') == 7
                                train_matlab = fullfile(train_dir, 'matlab');
                            else
                                train_matlab = train_dir;
                            end
                            addpath(genpath(train_matlab));
                            train_available = ~isempty(which('aps_linear')) && ...
                                              ~isempty(which('aps_weather_model')) && ...
                                              ~isempty(which('setparm_aps'));
                            if train_available
                                train_flag = 0;
                                updateOutput(app, ['TRAIN located and enabled for this run: ' train_matlab]);
                                persist_choice = uiconfirm(app.UIFigure, ...
                                    'TRAIN added for this session. Make it permanent so you do not have to locate it every time?', ...
                                    'TRAIN located', ...
                                    'Options', {'Make permanent', 'Just this session'}, ...
                                    'DefaultOption', 1, 'CancelOption', 2);
                                if strcmp(persist_choice, 'Make permanent')
                                    if savepath ~= 0
                                        up = userpath;
                                        if isempty(up)
                                            up = fullfile(char(java.lang.System.getProperty('user.home')), 'Documents', 'MATLAB');
                                        end
                                        if exist(up, 'dir') == 0
                                            mkdir(up);
                                        end
                                        fid_su = fopen(fullfile(up, 'startup.m'), 'a');
                                        if fid_su ~= -1
                                            fprintf(fid_su, '\naddpath(genpath(''%s''));\n', strrep(train_matlab, '\', '/'));
                                            fclose(fid_su);
                                            updateOutput(app, ['TRAIN path saved to ' fullfile(up, 'startup.m') ' - it will load at every MATLAB start.']);
                                        end
                                    else
                                        updateOutput(app, 'TRAIN path saved permanently (savepath).');
                                    end
                                end
                                break
                            else
                                retry_choice = uiconfirm(app.UIFigure, ...
                                    ['That folder does not contain TRAIN (aps_linear / aps_weather_model / setparm_aps were not found under it). ', ...
                                     'Choose a different folder, or skip the correction?'], ...
                                    'Not a TRAIN folder', ...
                                    'Options', {'Choose again...', 'Skip correction'}, ...
                                    'DefaultOption', 1, 'CancelOption', 2);
                                if strcmp(retry_choice, 'Skip correction')
                                    updateOutput(app, 'Proceeding without atmospheric correction.');
                                    break
                                end
                            end
                        end
                    else
                        updateOutput(app, 'Proceeding without atmospheric correction (user chose to skip).');
                    end
                end
                % end TRAIN interactive recovery
'''
END_SENTINEL_NEW = END_SENTINEL + RECOVERY


def edit(xml: str) -> str:
    if "% begin TRAIN interactive recovery" in xml:
        print("  locate-TRAIN dialog: already applied, skipping")
        return xml
    for label, before in (("intent snapshot", BEGIN_SENTINEL), ("end sentinel", END_SENTINEL)):
        n = xml.count(before)
        if n != 1:
            raise SystemExit(f"  locate-TRAIN dialog: anchor '{label}' found {n} times; expected 1. Aborting.")
    print("  locate-TRAIN dialog: applying")
    xml = xml.replace(BEGIN_SENTINEL, BEGIN_SENTINEL_NEW, 1)
    xml = xml.replace(END_SENTINEL, END_SENTINEL_NEW, 1)
    return xml


if __name__ == "__main__":
    edit_mlapp(MLAPP, edit)
    print("Done. Verify with: python tools/apply_train_locate_dialog.py  (should say 'already applied')")
