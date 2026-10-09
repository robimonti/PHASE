"""Extract PHASE_Preprocessing.mlapp into an editable beta engine class.

The stable MLAPP remains the source of truth while the beta is evaluated.  The
generated class preserves every callback and component so download, map,
update, import and processing code remain available.  Only class visibility,
repository-root resolution, the hidden engine window, the live process bridge
and the beta StaMPS handoff are adapted for the new controller.
"""

from pathlib import Path
import re
import zipfile


ROOT = Path(__file__).resolve().parents[1]
MLAPP = ROOT / "legacy" / "PHASE_Preprocessing.mlapp"
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


def _replace_checked(text: str, old: str, new: str, expected: int) -> str:
    count = text.count(old)
    if count != expected:
        raise RuntimeError(
            f"Expected {expected} processing anchor(s), found {count}: {old!r}"
        )
    return text.replace(old, new)


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
        + "        ExternalLogCallback              = []\n"
        + "        ExternalProgressCallback         = []\n"
        + "        ProjectRoot                      = phase_preprocessing_beta.projectRoot()\n"
        + "        InstallRoot                      = phase_preprocessing_beta.projectRoot()\n",
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
    # Legacy UI inventory callbacks are still reached by SaveButtonPushed.
    # Point their reads at the project data folder as well, even though the
    # beta controller presents its own inventory and downloader.
    code = _replace_checked(
        code, 'fullfile(appPath, "PHASE_Preprocessing", "slaves")',
        "fullfile(phase_preprocessing_beta.dataFolder(app.ProjectRoot), 'slaves')", 3,
    )
    code = _replace_checked(
        code, "fullfile(currentFolder, 'PHASE_Preprocessing', 'slaves')",
        "fullfile(phase_preprocessing_beta.dataFolder(app.ProjectRoot), 'slaves')", 6,
    )
    code = _replace_checked(
        code, "fullfile(currentFolder, 'PHASE_Preprocessing')",
        "phase_preprocessing_beta.dataFolder(app.ProjectRoot)", 5,
    )
    code = _replace_checked(
        code, 'fullfile(appPath, "downloadasf", "login_result.json")',
        'fullfile(app.ProjectRoot, "downloadasf", "login_result.json")', 3,
    )
    code = _replace_checked(
        code, 'fullfile(appPath, "downloadasf", "login_request.json")',
        'fullfile(app.ProjectRoot, "downloadasf", "login_request.json")', 2,
    )
    code = _replace_checked(
        code,
        "currentFolder = phase_preprocessing_beta.projectRoot();\n            cd(currentFolder);",
        "currentFolder = app.ProjectRoot;\n            cd(currentFolder);", 1,
    )
    code = _replace_checked(
        code, "function app = LegacyEngine\n",
        "function app = LegacyEngine(projectRoot, installRoot)\n"
        "            if nargin < 1 || isempty(projectRoot)\n"
        "                projectRoot = phase_preprocessing_beta.projectRoot();\n"
        "            end\n"
        "            if nargin < 2 || isempty(installRoot)\n"
        "                installRoot = phase_preprocessing_beta.projectRoot();\n"
        "            end\n"
        "            app.ProjectRoot = char(string(projectRoot));\n"
        "            app.InstallRoot = char(string(installRoot));\n", 1,
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

    # Keep code in the installation and all generated configuration, scripts,
    # SNAP products and StaMPS datasets in the selected project. Both sensor
    # branches use the same legacy anchors, so guard their occurrence counts.
    code = _replace_checked(code, "prep_folder = pwd;", "prep_folder = app.ProjectRoot;", 2)
    code = _replace_checked(code, "addpath(genpath(prep_folder));", "addpath(genpath(app.InstallRoot));", 2)
    code = _replace_checked(
        code,
        "load(strcat('.', par, 'PHASE_Preprocessing', par, 'input_preprocessing.mat'), 'python', ...",
        "load(fullfile(phase_preprocessing_beta.dataFolder(app.ProjectRoot), 'input_preprocessing.mat'), 'python', ...",
        2,
    )
    code = _replace_checked(
        code,
        "project_path_full = strcat(prep_folder, par, 'PHASE_Preprocessing');",
        "project_path_full = phase_preprocessing_beta.dataFolder(app.ProjectRoot);\n"
        "                        if ~isfolder(fullfile(project_path_full,'snap2stamps','bin'))\n"
        "                            mkdir(fullfile(project_path_full,'snap2stamps','bin'));\n"
        "                        end",
        2,
    )
    code = _replace_checked(
        code,
        "project_path_full, par, 'snap2stamps', par, 'graphs'",
        "app.InstallRoot, par, 'PHASE_Preprocessing', par, 'snap2stamps', par, 'graphs'",
        4,
    )
    code = _replace_checked(
        code, "project_parent_path_full = prep_folder;",
        "project_parent_path_full = phase_preprocessing_beta.stampsFolder(app.ProjectRoot);", 2,
    )
    code = _replace_checked(
        code, "filename = './PHASE_Preprocessing/input_preprocessing.mat';",
        "filename = fullfile(phase_preprocessing_beta.dataFolder(app.ProjectRoot), 'input_preprocessing.mat');", 2,
    )
    code = _replace_checked(
        code, "if exist('./PHASE_Preprocessing/input_preprocessing.mat', 'file') == 2",
        "if exist(fullfile(phase_preprocessing_beta.dataFolder(app.ProjectRoot), 'input_preprocessing.mat'), 'file') == 2", 1,
    )
    code = _replace_checked(
        code, "cfg = load('./PHASE_Preprocessing/input_preprocessing.mat');",
        "cfg = load(fullfile(phase_preprocessing_beta.dataFolder(app.ProjectRoot), 'input_preprocessing.mat'));", 1,
    )

    # A generated wrapper may live in the project, but it invokes the Python
    # source from the installation by absolute path and passes an absolute
    # per-project configuration path. This also handles spaces in paths.
    command_pattern = re.compile(
        r"(?P<name>step(?:_selector|1_master|[1-6]_slaves)) = "
        r"\('(?P<script>(?:SEN|CSK)_[A-Za-z0-9_]+\.py) "
        r"(?P<config>project(?:_master)?\.conf)'\);"
    )
    matches = list(command_pattern.finditer(code))
    if len(matches) != 16:
        raise RuntimeError(f"Expected 16 SNAP Python commands, found {len(matches)}")
    def command_replacement(match: re.Match[str]) -> str:
        return (f"{match['name']} = phase_preprocessing_beta.scriptCommand(python, "
                f"'{match['script']}', fullfile(project_path_full, 'snap2stamps', "
                f"'bin', '{match['config']}'), app.InstallRoot);")
    code = command_pattern.sub(command_replacement, code)
    code = code.replace("[python space step", "[step")
    code = code.replace("[dp python space step", "[dp step")
    code = code.replace("dp = ('::');", "dp = phase_preprocessing_beta.skipPrefix();")
    code = code.replace("[step_selector]", "step_selector")
    code = code.replace("[step1_master]", "step1_master")
    code = code.replace("step_selector_cmd = [step_selector];", "step_selector_cmd = step_selector;")
    code = code.replace("step_master_2 = [step1_master];", "step_master_2 = step1_master;")
    code = re.sub(r"= \[(step[1-6]_slaves)\];", r"= \1;", code)
    code = code.replace("[step5_slaves]", "step5_slaves")
    # ProcessBuilder executes .sh through bash, so chmod is unnecessary and
    # its unquoted path fails when a project name contains spaces.
    code = _replace_checked(code, "system(chmod);", "", 6)
    for handle, expected in (("f_snap2stamps_master", 4),
                             ("f_snap2stamps_slaves", 4),
                             ("f_average_intensity", 2)):
        code = _replace_checked(
            code,
            f"fprintf({handle},'#!/bin/bash \\n');",
            f"fprintf({handle},'#!/bin/bash \\nset -e\\n');",
            expected,
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

    # The copied MLAPP used to create the dataset directory as a side effect.
    # The text beta must create it before copying any diagnostic/configuration.
    stamps_folder_anchor = (
        "                        stamps_folder_full = "
        "fullfile(project_parent_path_full, stamps_folder);\n"
    )
    if code.count(stamps_folder_anchor) != 2:
        raise RuntimeError("Could not locate both StaMPS dataset folder assignments")
    code = code.replace(
        stamps_folder_anchor,
        stamps_folder_anchor
        + "                        if ~isfolder(stamps_folder_full)\n"
        + "                            mkdir(stamps_folder_full);\n"
        + "                            updateOutput(app, ['Created StaMPS dataset folder: ' stamps_folder_full]);\n"
        + "                        end\n",
    )

    # Open the text-based StaMPS beta from its canonical installation and pass
    # the new dataset folder explicitly.
    stable_launcher = "stamps_app_file = fullfile(stamps_app_full, 'PHASE_StaMPS.mlapp');"
    beta_launcher = "stamps_app_file = fullfile(project_path_full, 'PHASE_StaMPS_beta.m');"
    if code.count(stable_launcher) != 2:
        raise RuntimeError("Could not locate both stable StaMPS launch assignments")
    code = code.replace(stable_launcher, beta_launcher)
    code = _replace_checked(
        code, "stamps_app_file = fullfile(project_path_full, 'PHASE_StaMPS_beta.m');",
        "stamps_app_file = fullfile(app.InstallRoot, 'PHASE_Preprocessing', 'PHASE_StaMPS_beta.m');", 2,
    )
    code = _replace_checked(
        code, "stamps_diagnostic = fullfile(project_path_full, 'diagnose_PHASE_StaMPS.m');",
        "stamps_diagnostic = fullfile(app.InstallRoot, 'PHASE_Preprocessing', 'diagnose_PHASE_StaMPS.m');", 2,
    )
    if code.count("run(stamps_app_file);") != 2:
        raise RuntimeError("Could not locate both stable StaMPS run calls")
    code = code.replace(
        "run(stamps_app_file);",
        "phase_preprocessing_beta.launchStampsBeta(stamps_app_file, stamps_app_full);",
    )
    code = _replace_checked(
        code,
        """\
                                cd(stamps_app_full);
                                updateOutput(app, ['Changed MATLAB current directory to: ' stamps_app_full]);
""",
        """\
                                updateOutput(app, ['Launching the canonical PHASE_StaMPS_beta runtime for: ' stamps_app_full]);
""",
        2,
    )
    code = code.replace(
        "['Open it manually with: run(''' stamps_app_file ''')']",
        "['Open it manually with: PHASE_StaMPS_beta(''' stamps_app_full ''')']",
    )
    code = code.replace(
        "PHASE_StaMPS.mlapp not found in ' stamps_app_full",
        "PHASE_StaMPS_beta.m not found in ' project_path_full",
    )
    missing_stamps_config = """\
                            updateOutput(app, ['WARNING: input_StaMPS.mat was not found. ', ...
                                'Configure and save the StaMPS installation path before running.']);
"""
    bootstrap_stamps_config = """\
                            updateOutput(app, ['NOTICE: input_StaMPS.mat is not available yet. ', ...
                                'PHASE_StaMPS_beta will create it when the initial settings are saved.']);
"""
    code = _replace_checked(
        code,missing_stamps_config,bootstrap_stamps_config,2
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
                        updateOutput(app, ['Preprocessing completed. StaMPS dataset folder: ' stamps_app_full]);
                        choice = 'Open now';
                        if isfile(dst_input_mat)
                            updateOutput(app, ['Opening PHASE_StaMPS_beta with the existing configuration in: ' stamps_app_full]);
                        else
                            updateOutput(app, ['Opening PHASE_StaMPS_beta to create the initial configuration in: ' stamps_app_full]);
                        end
"""
        code = code[:start] + replacement + code[end:]
    code = code.replace(
        "% OPEN PHASE_StaMPS.mlapp (cross-platform).",
        "% OPEN PHASE_StaMPS_beta (cross-platform).",
    )

    # Replace platform terminals/asynchronous BAT execution with the beta's
    # silent process bridge. It runs the same generated scripts, but streams
    # their output and progress to the modern Run monitor.
    code = _replace_checked(
        code,
        "system(path_2_master);",
        "phase_preprocessing_beta.runCommandLive(app, path_2_master, ...\n"
        "                                    'Master selection and preparation', 3, 18, 1, 1);",
        2,
    )
    code = _replace_checked(
        code,
        "system(strjoin({chmod, path_2_master}, ';'));",
        "phase_preprocessing_beta.runCommandLive(app, path_2_master, ...\n"
        "                                    'Master selection and preparation', 3, 18, 1, 1);",
        4,
    )
    code = _replace_checked(
        code,
        "fullfile(project_path_full, '\\snap2stamps\\bin\\snap2stamps_slaves.bat &')",
        "fullfile(project_path_full, '\\snap2stamps\\bin\\snap2stamps_slaves.bat')",
        2,
    )
    code = _replace_checked(
        code,
        "                            path_2_slaves = [xterm space path_2_slaves];\n",
        "",
        2,
    )
    code = _replace_checked(
        code,
        "                            path_2_slaves = ['open -a Terminal ' path_2_slaves]; ",
        "                            % Executed without opening Terminal by runCommandLive.",
        1,
    )
    code = _replace_checked(
        code,
        "                            path_2_slaves = ['open -a Terminal ' path_2_slaves];",
        "                            % Executed without opening Terminal by runCommandLive.",
        1,
    )
    code = _replace_checked(
        code,
        "system(path_2_slaves); % execute the batch file",
        "phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...\n"
        "                                'Slave processing pipeline', 18, 92, first_step_num, 6);",
        2,
    )
    code = _replace_checked(
        code,
        "system(path_2_slaves, '-echo');",
        "phase_preprocessing_beta.runCommandLive(app, path_2_slaves, ...\n"
        "                                'Slave processing pipeline', 18, 92, first_step_num, 6);",
        4,
    )
    code = _replace_checked(
        code,
        "fullfile(project_path_full, '\\snap2stamps\\bin\\snap2stamps_update_average_intensity.bat &')",
        "fullfile(project_path_full, '\\snap2stamps\\bin\\snap2stamps_update_average_intensity.bat')",
        1,
    )
    code = _replace_checked(
        code,
        "                                path_2_average_intensity = [xterm space path_2_average_intensity];\n",
        "",
        1,
    )
    code = _replace_checked(
        code,
        "                                path_2_average_intensity = ['open -a Terminal ' path_2_average_intensity]; ",
        "                                % Executed without opening Terminal by runCommandLive.",
        1,
    )
    code = _replace_checked(
        code,
        "system(path_2_average_intensity);",
        "phase_preprocessing_beta.runCommandLive(app, path_2_average_intensity, ...\n"
        "                                    'Full-stack average intensity', 78, 92, 5, 5);",
        1,
    )
    code = _replace_checked(
        code,
        "system(path_2_average_intensity, '-echo');",
        "phase_preprocessing_beta.runCommandLive(app, path_2_average_intensity, ...\n"
        "                                    'Full-stack average intensity', 78, 92, 5, 5);",
        2,
    )
    code = _replace_checked(
        code,
        "                        updateOutput(app, 'An error occurred during script execution. Please check each step log file!');",
        "                        if strcmp(ME.identifier,'PHASE:ProcessingStopped')\n"
        "                            updateOutput(app, 'Preprocessing was force-stopped by the user.');\n"
        "                        else\n"
        "                            updateOutput(app, 'An error occurred during script execution. Please check each step log file!');\n"
        "                        end",
        2,
    )

    code = _replace_last(
        code,
        "app.UIFigure.Visible = 'on';",
        "app.UIFigure.Visible = 'off';",
    )

    # The hub owns navigation between modules. Keep the exported dataset, but
    # never open a second StaMPS window after preprocessing finishes.
    launch_start = "                        % OPEN PHASE_StaMPS_beta (cross-platform)."
    launch_end = "                        %% ----------------------------------------------------"
    if code.count(launch_start) != 2:
        raise RuntimeError("Could not locate both StaMPS auto-launch blocks")
    for _ in range(2):
        start = code.index(launch_start)
        end = code.index(launch_end, start)
        replacement = """\
                        stamps_app_full = fullfile(project_parent_path_full, stamps_folder);
                        updateOutput(app, ['Preprocessing completed. StaMPS dataset folder: ' stamps_app_full]);
                        updateOutput(app, 'Open the StaMPS PSI section in the PHASE hub when ready.');

"""
        code = code[:start] + replacement + code[end:]

    header = """\
% GENERATED FROM PHASE_Preprocessing.mlapp.
% Run tools/extract_preprocessing_beta_engine.py after changing the stable app.
% The complete formerly embedded App Designer source follows as editable text.

"""
    # App Designer XML contains pervasive indentation-only/trailing spaces.
    # Normalise them so the editable generated backend remains diff-clean.
    code = "\n".join(line.rstrip() for line in code.splitlines())
    return header + code.rstrip() + "\n"


def main() -> None:
    with zipfile.ZipFile(MLAPP) as archive:
        xml = archive.read("matlab/document.xml").decode("utf-8")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(extract(xml), encoding="utf-8", newline="\n")
    print(f"Generated {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
