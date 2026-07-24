"""Extract PHASE_Preprocessing.mlapp into an editable beta engine class.

The stable MLAPP remains the source of truth while the beta is evaluated.  The
generated class preserves every callback and component so download, map,
update, import and processing code remain available.  Only class visibility,
repository-root resolution, the hidden engine window and an external log hook
are adapted for the new controller.
"""

from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]
MLAPP = ROOT / "PHASE_Preprocessing.mlapp"
OUTPUT = (
    ROOT
    / "PHASE_Preprocessing"
    / "+phase_preprocessing_beta"
    / "LegacyEngine.m"
)


def _replace_last(text: str, old: str, new: str) -> str:
    index = text.rfind(old)
    if index < 0:
        raise RuntimeError(f"Could not locate final engine anchor: {old!r}")
    return text[:index] + new + text[index + len(old):]


def extract(xml: str) -> str:
    cdata_start = xml.index("<![CDATA[") + len("<![CDATA[")
    cdata_end = xml.rindex("]]>")
    code = xml[cdata_start:cdata_end].replace("\r\n", "\n")

    code = code.replace(
        "classdef PHASE_Preprocessing < matlab.apps.AppBase",
        "classdef LegacyEngine < matlab.apps.AppBase",
        1,
    )
    code = code.replace(
        "function app = PHASE_Preprocessing",
        "function app = LegacyEngine",
        1,
    )
    code = code.replace("properties (Access = private)", "properties (Access = public)")
    code = code.replace("methods (Access = private)", "methods (Access = public)")

    component_anchor = """\
        ContextMenu                     matlab.ui.container.ContextMenu
        Menu                            matlab.ui.container.Menu
        Menu2                           matlab.ui.container.Menu
"""
    if code.count(component_anchor) != 1:
        raise RuntimeError("Could not locate the public component property block")
    code = code.replace(
        component_anchor,
        component_anchor
        + "        ExternalLogCallback              = []\n",
        1,
    )

    # The MLAPP lived in the repository root; the generated class lives two
    # folders deeper.  Keep every path calculation semantically identical.
    code = code.replace(
        "fileparts(mfilename('fullpath'))",
        "phase_preprocessing_beta.projectRoot()",
    )
    code = code.replace(
        'fileparts(mfilename("fullpath"))',
        "phase_preprocessing_beta.projectRoot()",
    )

    output_anchor = """\
            end
            
        end

        % Downloader parameters
"""
    output_hook = """\
            end

            if ~isempty(app.ExternalLogCallback)
                try
                    app.ExternalLogCallback(message);
                catch callbackError
                    warning('PHASE:PreprocessingBetaLogCallback', ...
                        'Could not forward engine log: %s', callbackError.message);
                end
            end
            
        end

        % Downloader parameters
"""
    if code.count(output_anchor) != 1:
        raise RuntimeError("Could not locate updateOutput callback boundary")
    code = code.replace(output_anchor, output_hook, 1)

    code = code.replace(
        "app.UIFigure.Name = 'MATLAB App';",
        "app.UIFigure.Name = 'PHASE · Preprocessing advanced workspace';",
        1,
    )
    code = code.replace(
        "% Set the working directory to the folder where the .mlapp file is located",
        "% Set the working directory to the PHASE project root",
        1,
    )

    # The stable app still copies PHASE_StaMPS.mlapp into every ASC/DES folder.
    # The beta must have no runtime MLAPP dependency: remove those two copy
    # blocks and retain only the optional diagnostic helper beside the data.
    legacy_copy_start = "                        dest_mlapp = fullfile(stamps_folder_full, 'PHASE_StaMPS.mlapp');"
    legacy_copy_end = "                        % EVENTUAL REMOVAL OF THE SLAVES IMAGES TO SAVE DISK SPACE"
    if code.count(legacy_copy_start) != 2 or code.count(legacy_copy_end) != 2:
        raise RuntimeError("Could not locate both stable StaMPS copy blocks")
    for _ in range(2):
        start = code.index(legacy_copy_start)
        end = code.index(legacy_copy_end, start)
        replacement = """\
                        % The beta launches the editable StaMPS module from its
                        % canonical installation. No PHASE_StaMPS.mlapp is copied
                        % or required at runtime.
                        stamps_diagnostic = fullfile(project_path_full, 'diagnose_PHASE_StaMPS.m');
                        if isfile(stamps_diagnostic)
                            copyfile(stamps_diagnostic, stamps_folder_full, 'f');
                        end
                        
"""
        code = code[:start] + replacement + code[end:]

    # Open the text-based StaMPS beta from its canonical installation and pass
    # the new dataset folder explicitly.
    stable_launcher = "stamps_app_file = fullfile(stamps_app_full, 'PHASE_StaMPS.mlapp');"
    beta_launcher = "stamps_app_file = fullfile(project_path_full, 'PHASE_StaMPS_beta.m');"
    if code.count(stable_launcher) != 2:
        raise RuntimeError("Could not locate both stable StaMPS launch assignments")
    code = code.replace(stable_launcher, beta_launcher)
    if code.count("run(stamps_app_file);") != 2:
        raise RuntimeError("Could not locate both stable StaMPS run calls")
    code = code.replace(
        "run(stamps_app_file);",
        "PHASE_StaMPS_beta(stamps_app_full);",
    )
    code = code.replace(
        "['Open it manually with: run(''' stamps_app_file ''')']",
        "['Open it manually with: PHASE_StaMPS_beta(''' stamps_app_full ''')']",
    )
    code = code.replace(
        "PHASE_StaMPS.mlapp not found in ' stamps_app_full",
        "PHASE_StaMPS_beta.m not found in ' project_path_full",
    )
    prompt_start = """\
                        try
                            choice = uiconfirm(app.UIFigure, ...
"""
    prompt_end = "                        if strcmp(choice, 'Open now')"
    if code.count(prompt_start) != 2:
        raise RuntimeError("Could not locate both stable StaMPS confirmation dialogs")
    for _ in range(2):
        start = code.index(prompt_start)
        end = code.index(prompt_end, start)
        replacement = """\
                        choice = 'Open now';
                        updateOutput(app, ['Preprocessing completed. Opening PHASE_StaMPS_beta in: ' stamps_app_full]);
"""
        code = code[:start] + replacement + code[end:]
    code = code.replace(
        "% OPEN PHASE_StaMPS.mlapp (cross-platform).",
        "% OPEN PHASE_StaMPS_beta (cross-platform).",
    )
    code = _replace_last(
        code,
        "app.UIFigure.Visible = 'on';",
        "app.UIFigure.Visible = 'off';",
    )

    header = """\
% GENERATED FROM PHASE_Preprocessing.mlapp.
% Run tools/extract_preprocessing_beta_engine.py after changing the stable app.
% The complete formerly embedded App Designer source follows as editable text.

"""
    return header + code.rstrip() + "\n"


def main() -> None:
    with zipfile.ZipFile(MLAPP) as archive:
        xml = archive.read("matlab/document.xml").decode("utf-8")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(extract(xml), encoding="utf-8", newline="\n")
    print(f"Generated {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
