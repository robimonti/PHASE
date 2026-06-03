"""Add a 'TS Points' tab + 'Open TS Points picker...' button to
PHASE_Preprocessing/PHASE_StaMPS.mlapp. Idempotent: re-running detects
the new component already exists and exits with a no-op.

The .mlapp is an OOXML zip; this script edits matlab/document.xml in
place via tools/mlapp_roundtrip.py and never touches appdesigner/appModel.mat.
App Designer will still open the file but treats the new components as
"unknown" until you re-save them through the GUI; the runtime always
runs the patched code, so functionally everything works.
"""
from pathlib import Path
import sys

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))
from mlapp_roundtrip import edit_mlapp  # noqa: E402

MLAPP = THIS_DIR.parent / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"

PROPERTY_DECLS = """\
        TSPointsTab                     matlab.ui.container.Tab
        OpenTSPickerButton              matlab.ui.control.Button
"""

STARTUP_ADDPATH = """\

            % Wire StaMPS matlab + matlab_compat onto the path so the
            % TS Points picker (ts_export_picker / ts_export_batch /
            % llh2local) resolves regardless of the user's startup.m.
            %
            % NOTE: PathEditField_2 = StaMPS installation folder (label
            % StaMPSinstallationfolderLabel), PathEditField_3 = Project
            % folder (label ProjectfolderLabel). These names are
            % inherited from the existing PHASE_StaMPS layout. If the
            % .mlapp tabs are reordered, also update the references in
            % OpenTSPickerButtonPushed below.
            stamps_root = app.PathEditField_2.Value;
            if ~isempty(stamps_root) && isfolder(stamps_root)
                addpath(fullfile(stamps_root, 'matlab'));
                addpath(fullfile(stamps_root, 'matlab_compat'));
            end
"""

TAB_CREATION = """\

            % Create TSPointsTab
            app.TSPointsTab = uitab(app.TabGroup);
            app.TSPointsTab.Title = 'TS Points';
            app.TSPointsTab.BackgroundColor = [1 1 1];

            % Welcome content (removed and replaced by the embedded
            % picker once OpenTSPickerButton is pushed).
            tsLabel = uilabel(app.TSPointsTab);
            tsLabel.FontName = 'Manrope';
            tsLabel.FontWeight = 'bold';
            tsLabel.FontSize = 14;
            tsLabel.Position = [40 572 1100 22];
            tsLabel.Text = ['Pick query points (manually or by clicking on the velocity map) ' ...
                            'and export their time series as CSV.'];

            tsLabel2 = uilabel(app.TSPointsTab);
            tsLabel2.FontName = 'Manrope';
            tsLabel2.Position = [40 540 1100 22];
            tsLabel2.Text = ['Requires stamps(7) + ps_plot(''v-do'',''ts'',1,...) to have produced ' ...
                             'ps_plot_ts_v-do.mat in the project folder set on the Preparation tab.'];

            tsLabel3 = uilabel(app.TSPointsTab);
            tsLabel3.FontName = 'Manrope';
            tsLabel3.FontAngle = 'italic';
            tsLabel3.Position = [40 510 1100 22];
            tsLabel3.Text = ['If <project folder>/aoi_points.csv already exists, the picker auto-loads it on start.'];

            % Create OpenTSPickerButton — clicking it embeds the picker
            % into THIS tab (no popup window).
            app.OpenTSPickerButton = uibutton(app.TSPointsTab, 'push');
            app.OpenTSPickerButton.ButtonPushedFcn = createCallbackFcn(app, @OpenTSPickerButtonPushed, true);
            app.OpenTSPickerButton.BackgroundColor = [0.4 0.7 0.4];
            app.OpenTSPickerButton.FontName = 'Manrope';
            app.OpenTSPickerButton.FontSize = 14;
            app.OpenTSPickerButton.Position = [475 460 280 44];
            app.OpenTSPickerButton.Text = 'Load TS picker into this tab';
"""

CALLBACK_METHOD = """\

        % WARNING: TS Points tab + OpenTSPickerButton + TSPointsTab
        % property + this callback are added by tools/add_tspoints_tab.py
        % which patches matlab/document.xml only, NOT
        % appdesigner/appModel.mat. Re-saving this app from App Designer
        % may strip them. Re-run the patcher to restore.
        % See docs/TS_export_picker_integration.md.
        %
        % Button pushed function: OpenTSPickerButton
        % Embeds the time-series picker INSIDE app.TSPointsTab — the
        % welcome labels and the button itself are wiped (the picker
        % calls delete(parent.Children)).
        function OpenTSPickerButtonPushed(app, event)
            workdir = app.PathEditField_3.Value;
            if isempty(workdir) || ~isfolder(workdir)
                uialert(app.UIFigure, ...
                    'Set the project folder in the Preparation tab first.', ...
                    'Workdir missing');
                return;
            end
            try
                ts_export_picker(workdir, app.TSPointsTab);
            catch ME
                uialert(app.UIFigure, ME.message, 'Picker failed');
            end
        end
"""


def patch(xml: str) -> str:
    if "OpenTSPickerButton" in xml:
        print("already patched, no-op")
        return xml

    # 1. Property declarations: insert right after the public properties opening line.
    needle = "    properties (Access = public)\n        UIFigure"
    replacement = (
        "    properties (Access = public)\n"
        + PROPERTY_DECLS
        + "        UIFigure"
    )
    xml = xml.replace(needle, replacement, 1)
    if "TSPointsTab " not in xml:
        raise RuntimeError("property insertion failed (needle not found)")

    # 2. startupFcn: append addpath right after `cd(currentFolder);`.
    cd_needle = "            cd(currentFolder);\n"
    if cd_needle not in xml:
        raise RuntimeError("startupFcn cd line not found")
    xml = xml.replace(cd_needle, cd_needle + STARTUP_ADDPATH, 1)

    # 3. Tab + button creation: insert before "Show the figure after all
    #    components are created" (last block of createComponents).
    show_needle = "            % Show the figure after all components are created"
    if show_needle not in xml:
        raise RuntimeError("createComponents end marker not found")
    xml = xml.replace(show_needle, TAB_CREATION + "\n" + show_needle, 1)

    # 4. Callback method: the callbacks `methods (Access = private)`
    #    block closes with `    end` immediately before the comment
    #    `    % App creation and deletion`. Inject the new callback
    #    just before that closing `end`, inside the same block.
    cb_close = "    end\n\n    % App creation and deletion"
    if cb_close not in xml:
        raise RuntimeError("callbacks block close marker not found")
    xml = xml.replace(cb_close, CALLBACK_METHOD + cb_close, 1)

    return xml


if __name__ == "__main__":
    if not MLAPP.exists():
        sys.exit(f"missing: {MLAPP}")
    edit_mlapp(MLAPP, patch)
    print(f"patched: {MLAPP}")
