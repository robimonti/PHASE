"""Keep the canonical StaMPS app and copy only the current version to datasets.

PHASE_Preprocessing historically moved PHASE_StaMPS.mlapp out of the installed
PHASE_Preprocessing folder into the first ASC_/DES_ dataset.  That broke the
installer shortcut; subsequent datasets then copied whichever app happened to
exist in an older dataset, allowing stale Windows code to propagate forever.

This migration updates both Sentinel-1 and COSMO-SkyMed branches in the binary
App Designer archive. The canonical app remains installed and is force-copied
to the current dataset on every completed preprocessing run. A dataset's
existing input_StaMPS.mat is deliberately preserved because it may contain
user-tuned processing parameters; the canonical MAT is copied only when the
dataset does not have one yet.
"""

import re
from pathlib import Path

from mlapp_roundtrip import edit_mlapp


REPO_ROOT = Path(__file__).resolve().parents[1]
MLAPP = REPO_ROOT / "PHASE_Preprocessing.mlapp"


HANDOFF_PATTERN = re.compile(
    r"                        project_parent_path_full = prep_folder;\n"
    r"                        exist_StaMPS_A = exist\(strcat\(project_path_full, par, 'PHASE_StaMPS\.mlapp'\), \"file\"\);\n"
    r".*?"
    r"                        end\n"
    r"                        \n"
    r"(?=                        % EVENTUAL REMOVAL OF THE SLAVES IMAGES TO SAVE DISK SPACE)",
    re.DOTALL,
)


HANDOFF_REPLACEMENT = """\
                        project_parent_path_full = prep_folder;
                        stamps_folder_full = fullfile(project_parent_path_full, stamps_folder);
                        dest_mlapp = fullfile(stamps_folder_full, 'PHASE_StaMPS.mlapp');
                        canonical_stamps_app = fullfile(project_path_full, 'PHASE_StaMPS.mlapp');

                        if ~isfolder(stamps_folder_full)
                            mkdir(stamps_folder_full);
                            updateOutput(app, ['Created new StaMPS folder: ' stamps_folder]);
                        else
                            updateOutput(app, 'The StaMPS processing folder already exists - its launcher will be refreshed and its saved parameters preserved.');
                        end

                        % The installed app is the single source of truth.  Keep it in
                        % PHASE_Preprocessing so the installer never leaves a broken
                        % shortcut, and overwrite dataset copies so an old ASC_/DES_
                        % folder can never propagate stale Windows/source code.
                        if isfile(canonical_stamps_app)
                            copyfile(canonical_stamps_app, dest_mlapp, 'f');
                            updateOutput(app, ['Copied current PHASE_StaMPS.mlapp to: ' stamps_folder_full]);
                            % Keep the one-command runtime diagnostic beside each
                            % dataset app so Windows failures can be reported without
                            % asking users to transcribe several Command Window calls.
                            stamps_diagnostic = fullfile(project_path_full, 'diagnose_PHASE_StaMPS.m');
                            if isfile(stamps_diagnostic)
                                copyfile(stamps_diagnostic, stamps_folder_full, 'f');
                            end
                        else
                            error('PHASE_Preprocessing:StaMPSAppMissing', ...
                                ['The canonical PHASE_StaMPS.mlapp is missing from %s. ', ...
                                 'Repair/update the PHASE installation; a copy from an older ', ...
                                 'ASC_/DES_ folder will not be used because it may be stale.'], ...
                                 project_path_full);
                        end
                        
"""


INPUT_PATTERN = re.compile(
    r"                        % Bring along input_StaMPS\.mat .*?"
    r"(?=                        % OPEN PHASE_StaMPS\.mlapp \(cross-platform\)\.)",
    re.DOTALL,
)


INPUT_REPLACEMENT = """\
                        % Seed a new dataset with the installed configuration, but never
                        % overwrite an existing input_StaMPS.mat: it may contain parameters
                        % the user tuned specifically for this processing run.
                        src_input_mat = fullfile(project_path_full, 'input_StaMPS.mat');
                        dst_input_mat = fullfile(stamps_folder_full, 'input_StaMPS.mat');
                        if isfile(dst_input_mat)
                            updateOutput(app, 'Preserved the existing dataset input_StaMPS.mat.');
                        elseif isfile(src_input_mat)
                            copyfile(src_input_mat, dst_input_mat);
                            updateOutput(app, 'Copied input_StaMPS.mat to the new StaMPS processing folder.');
                        else
                            updateOutput(app, ['WARNING: input_StaMPS.mat was not found. ', ...
                                'Configure and save the StaMPS installation path before running.']);
                        end

"""


def patch(xml: str) -> str:
    marker = "The installed app is the single source of truth."
    if xml.count(marker) == 2:
        changed = False
        unsafe_input = """\
                        % Refresh input_StaMPS.mat together with the app.  Keeping an
                        % older file in an existing dataset can point MATLAB at a stale
                        % StaMPS/TRAIN installation after an installer update.
                        src_input_mat = fullfile(project_path_full, 'input_StaMPS.mat');
                        dst_input_mat = fullfile(stamps_folder_full, 'input_StaMPS.mat');
                        if isfile(src_input_mat)
                            copyfile(src_input_mat, dst_input_mat, 'f');
                            updateOutput(app, 'Refreshed input_StaMPS.mat in the StaMPS processing folder.');
                        elseif isfile(dst_input_mat)
                            updateOutput(app, ['WARNING: canonical input_StaMPS.mat is missing; ', ...
                                'the existing dataset copy was left unchanged. Verify the StaMPS installation path before running.']);
                        else
                            updateOutput(app, ['WARNING: input_StaMPS.mat was not found. ', ...
                                'Configure and save the StaMPS installation path before running.']);
                        end

"""
        if xml.count(unsafe_input) == 2:
            xml = xml.replace(unsafe_input, INPUT_REPLACEMENT)
            xml = xml.replace(
                "The StaMPS processing folder already exists - its launcher and parameters will be refreshed.",
                "The StaMPS processing folder already exists - its launcher will be refreshed and its saved parameters preserved.",
            )
            changed = True

        diagnostic_marker = "Keep the one-command runtime diagnostic beside each"
        if xml.count(diagnostic_marker) == 2:
            if changed:
                print("Updated both handoffs to preserve dataset input_StaMPS.mat files.")
            else:
                print("PHASE_Preprocessing.mlapp already contains the StaMPS handoff fix.")
            return xml
        copy_anchor = (
            "                            updateOutput(app, ['Copied current "
            "PHASE_StaMPS.mlapp to: ' stamps_folder_full]);\n"
        )
        diagnostic_block = copy_anchor + (
            "                            % Keep the one-command runtime diagnostic beside each\n"
            "                            % dataset app so Windows failures can be reported without\n"
            "                            % asking users to transcribe several Command Window calls.\n"
            "                            stamps_diagnostic = fullfile(project_path_full, 'diagnose_PHASE_StaMPS.m');\n"
            "                            if isfile(stamps_diagnostic)\n"
            "                                copyfile(stamps_diagnostic, stamps_folder_full, 'f');\n"
            "                            end\n"
        )
        if xml.count(copy_anchor) != 2:
            raise RuntimeError("Could not upgrade both existing StaMPS handoff blocks")
        print("Added the StaMPS diagnostic helper to both dataset handoffs.")
        return xml.replace(copy_anchor, diagnostic_block)

    handoff_matches = list(HANDOFF_PATTERN.finditer(xml))
    input_matches = list(INPUT_PATTERN.finditer(xml))
    if len(handoff_matches) != 2:
        raise RuntimeError(
            f"Expected Sentinel-1 and COSMO handoff blocks (2), found {len(handoff_matches)}"
        )
    if len(input_matches) != 2:
        raise RuntimeError(
            f"Expected Sentinel-1 and COSMO input blocks (2), found {len(input_matches)}"
        )

    xml = HANDOFF_PATTERN.sub(HANDOFF_REPLACEMENT, xml)
    return INPUT_PATTERN.sub(INPUT_REPLACEMENT, xml)


if __name__ == "__main__":
    if not MLAPP.is_file():
        raise SystemExit(f"Expected mlapp at {MLAPP}")
    edit_mlapp(MLAPP, patch)
    print(f"Patched {MLAPP.relative_to(REPO_ROOT)}")
