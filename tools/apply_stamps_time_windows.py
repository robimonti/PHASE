"""Split the three StaMPS temporal windows from the dataset time span.

PHASE historically passed the acquisition span to weed_time_win,
unwrap_time_win and scn_time_win. They are independent processing parameters.
Add one GUI control for each, persist them in input_StaMPS.mat, and retain old
MAT-file behavior by falling back to time_span when a field is absent.
"""

from pathlib import Path

from mlapp_roundtrip import edit_mlapp


REPO_ROOT = Path(__file__).resolve().parents[1]
MLAPP = REPO_ROOT / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"


def replace_once(xml: str, old: str, new: str, description: str) -> str:
    if xml.count(old) != 1:
        raise RuntimeError(f"{description} anchor was not found exactly once")
    return xml.replace(old, new, 1)


def patch(xml: str) -> str:
    marker = "% begin independent StaMPS temporal windows"
    app_designer_warning = """\
        % WARNING: weed_time_win, unwrap_time_win and scn_time_win GUI
        % controls are added by tools/apply_stamps_time_windows.py, which
        % patches matlab/document.xml only, NOT appdesigner/appModel.mat.
        % Re-saving this app from App Designer may strip them. Re-run the
        % patcher to restore them.

"""
    warning_anchor = "        % WARNING: TS Points tab + OpenTSPickerButton + TSPointsTab\n"
    if marker in xml:
        # Do not populate GUI controls from Start: doing so could silently
        # erase an unsaved edit. The runtime variables still fall back to
        # time_span, while the existing mismatch guard requests Load/Save.
        repaired = xml
        for name in ("weed_time_win", "unwrap_time_win", "scn_time_win"):
            repaired = repaired.replace(
                f"                    app.{name}EditField.Value = {name};\n", ""
            )
        repaired = repaired.replace(
            "% older value and could unexpectedly launch GACOS or seek tca2.",
            "% older value, launch GACOS/seek tca2, or apply a wrong temporal window.",
        )
        if app_designer_warning not in repaired:
            repaired = replace_once(
                repaired,
                warning_anchor,
                app_designer_warning + warning_anchor,
                "App Designer warning",
            )
        if repaired != xml:
            print("Repaired legacy time-window loading to preserve unsaved-edit detection.")
        else:
            print("PHASE_StaMPS.mlapp already contains independent temporal windows.")
        return repaired

    # Public App Designer component properties.
    xml = replace_once(
        xml,
        "        Step4PSweedingLabel             matlab.ui.control.Label\n",
        """\
        Step4PSweedingLabel             matlab.ui.control.Label
        weed_time_winDescriptionLabel   matlab.ui.control.Label
        weed_time_winEditField          matlab.ui.control.NumericEditField
        weed_time_winEditFieldLabel     matlab.ui.control.Label
""",
        "weed-time component property",
    )
    xml = replace_once(
        xml,
        "        Step6PhaseunwrappingLabel       matlab.ui.control.Label\n",
        """\
        Step6PhaseunwrappingLabel       matlab.ui.control.Label
        unwrap_time_winDescriptionLabel  matlab.ui.control.Label
        unwrap_time_winEditField        matlab.ui.control.NumericEditField
        unwrap_time_winEditFieldLabel   matlab.ui.control.Label
""",
        "unwrap-time component property",
    )
    xml = replace_once(
        xml,
        "        Step8AtmosphericfilteringLabel  matlab.ui.control.Label\n",
        """\
        Step8AtmosphericfilteringLabel  matlab.ui.control.Label
        scn_time_winDescriptionLabel    matlab.ui.control.Label
        scn_time_winEditField           matlab.ui.control.NumericEditField
        scn_time_winEditFieldLabel      matlab.ui.control.Label
""",
        "SCN-time component property",
    )

    # Independent defaults match ps_parms_default.m in StaMPS.
    xml = replace_once(
        xml,
        "        time_span = 0;\t\t\n",
        """\
        time_span = 0;		
        weed_time_win = 730;
        unwrap_time_win = 730;
        scn_time_win = 365;
""",
        "private temporal defaults",
    )

    # Save callback and MAT persistence.
    xml = replace_once(
        xml,
        "            time_span = app.DaysEditField.Value;\t\t\n",
        """\
            time_span = app.DaysEditField.Value;		
            weed_time_win = app.weed_time_winEditField.Value;
            unwrap_time_win = app.unwrap_time_winEditField.Value;
            scn_time_win = app.scn_time_winEditField.Value;
""",
        "Save callback temporal controls",
    )
    mat_list_anchor = "'master_date', 'export_name', 'time_span', ..."
    if xml.count(mat_list_anchor) != 2:
        raise RuntimeError("Expected Save and Start MAT variable-list anchors")
    xml = xml.replace(
        mat_list_anchor,
        "'master_date', 'export_name', 'time_span', 'weed_time_win', 'unwrap_time_win', 'scn_time_win', ...",
    )

    # Old MAT files do not have these fields. Preserve their exact historical
    # behavior, surface it in the UI, and let the next explicit Save migrate.
    load_end = """\
                    'scla_method', 'scla_drop_index', 'scn_wavelength', 'scn_kriging_flag', 'ref_centre_lonlat', ...
                    'ref_radius', 'ref_velocity', 'plot_s', 'ref_centre_lonlat_w', 'ref_radius_w', 'ph_output');

                % begin saved/visible tropospheric-configuration guard
"""
    load_end_new = """\
                    'scla_method', 'scla_drop_index', 'scn_wavelength', 'scn_kriging_flag', 'ref_centre_lonlat', ...
                    'ref_radius', 'ref_velocity', 'plot_s', 'ref_centre_lonlat_w', 'ref_radius_w', 'ph_output');

                % begin independent StaMPS temporal windows
                legacy_time_window_fields = {};
                if exist('weed_time_win', 'var') ~= 1
                    weed_time_win = time_span;
                    legacy_time_window_fields{end + 1} = 'weed_time_win';
                end
                if exist('unwrap_time_win', 'var') ~= 1
                    unwrap_time_win = time_span;
                    legacy_time_window_fields{end + 1} = 'unwrap_time_win';
                end
                if exist('scn_time_win', 'var') ~= 1
                    scn_time_win = time_span;
                    legacy_time_window_fields{end + 1} = 'scn_time_win';
                end
                if ~isempty(legacy_time_window_fields)
                    updateOutput(app, ['Legacy input_StaMPS.mat: missing ' ...
                        strjoin(legacy_time_window_fields, ', ') ...
                        '; using time_span=' num2str(time_span) ...
                        ' days for backward compatibility. Press Save to migrate the MAT file.']);
                end
                % end independent StaMPS temporal windows

                % begin saved/visible tropospheric-configuration guard
"""
    xml = replace_once(xml, load_end, load_end_new, "Start MAT-load end")

    visible_anchor = """\
                visible_train_flag = double(~app.TRAINatmosphericcorrectionCheckBox.Value);
                tropo_config_mismatch = ...
"""
    visible_new = """\
                visible_train_flag = double(~app.TRAINatmosphericcorrectionCheckBox.Value);
                visible_weed_time_win = app.weed_time_winEditField.Value;
                visible_unwrap_time_win = app.unwrap_time_winEditField.Value;
                visible_scn_time_win = app.scn_time_winEditField.Value;
                tropo_config_mismatch = ...
"""
    xml = replace_once(xml, visible_anchor, visible_new, "visible configuration values")

    condition_anchor = """\
                    train_flag ~= visible_train_flag;
                if tropo_config_mismatch
"""
    condition_new = """\
                    train_flag ~= visible_train_flag;
                time_window_config_mismatch = ...
                    ~isequaln(weed_time_win, visible_weed_time_win) || ...
                    ~isequaln(unwrap_time_win, visible_unwrap_time_win) || ...
                    ~isequaln(scn_time_win, visible_scn_time_win);
                if tropo_config_mismatch || time_window_config_mismatch
"""
    xml = replace_once(xml, condition_anchor, condition_new, "configuration mismatch guard")

    error_anchor = """\
                        ['The visible TRAIN/tropospheric controls do not match input_StaMPS.mat. ', ...
                         'Processing was not started and the MAT file was not changed. ', ...
                         'Press Load, make the intended changes, press Save, then press Start. ', ...
                         'Saved: TRAIN=%d, subtr_tropo=%s, tropo_method=%s. ', ...
                         'Visible: TRAIN=%d, subtr_tropo=%s, tropo_method=%s.'], ...
                         train_flag == 0, loaded_subtr_tropo, loaded_tropo_method, ...
                         visible_train_flag == 0, visible_subtr_tropo, visible_tropo_method);
"""
    error_new = """\
                        ['The visible TRAIN/tropospheric/temporal-window controls do not match input_StaMPS.mat. ', ...
                         'Processing was not started and the MAT file was not changed. ', ...
                         'Press Load, make the intended changes, press Save, then press Start. ', ...
                         'Saved: TRAIN=%d, subtr_tropo=%s, tropo_method=%s, time windows=%g/%g/%g days. ', ...
                         'Visible: TRAIN=%d, subtr_tropo=%s, tropo_method=%s, time windows=%g/%g/%g days.'], ...
                         train_flag == 0, loaded_subtr_tropo, loaded_tropo_method, ...
                         weed_time_win, unwrap_time_win, scn_time_win, ...
                         visible_train_flag == 0, visible_subtr_tropo, visible_tropo_method, ...
                         visible_weed_time_win, visible_unwrap_time_win, visible_scn_time_win);
"""
    xml = replace_once(xml, error_anchor, error_new, "unsaved configuration error")

    status_anchor = """\
                updateOutput(app, sprintf( ...
                    'Loaded configuration: TRAIN=%d, subtr_tropo=%s, tropo_method=%s', ...
                    train_flag == 0, loaded_subtr_tropo, loaded_tropo_method));
                clear loaded_subtr_tropo visible_subtr_tropo loaded_tropo_method
                clear visible_tropo_method visible_train_flag tropo_config_mismatch
"""
    status_new = """\
                updateOutput(app, sprintf( ...
                    ['Loaded configuration: TRAIN=%d, subtr_tropo=%s, tropo_method=%s, ', ...
                     'weed/unwrap/scn time windows=%g/%g/%g days'], ...
                    train_flag == 0, loaded_subtr_tropo, loaded_tropo_method, ...
                    weed_time_win, unwrap_time_win, scn_time_win));
                clear loaded_subtr_tropo visible_subtr_tropo loaded_tropo_method
                clear visible_tropo_method visible_train_flag tropo_config_mismatch
                clear visible_weed_time_win visible_unwrap_time_win visible_scn_time_win
                clear time_window_config_mismatch legacy_time_window_fields
"""
    xml = replace_once(xml, status_anchor, status_new, "loaded configuration status")

    # Runtime parameter mapping: no temporal processing setting may use the
    # acquisition span after this migration.
    runtime_replacements = {
        "setparm('weed_time_win', time_span);": "setparm('weed_time_win', weed_time_win);",
        "setparm('unwrap_time_win', time_span);": "setparm('unwrap_time_win', unwrap_time_win);",
        "setparm('scn_time_win', time_span);": "setparm('scn_time_win', scn_time_win);",
    }
    for old, new in runtime_replacements.items():
        xml = replace_once(xml, old, new, old)

    # Load callback, including a fallback for every legacy MAT file.
    load_ui_anchor = "            app.DaysEditField.Value = data.time_span;\t\t\n"
    load_ui_new = """\
            app.DaysEditField.Value = data.time_span;		
                if isfield(data, 'weed_time_win')
                    app.weed_time_winEditField.Value = data.weed_time_win;
                else
                    app.weed_time_winEditField.Value = data.time_span;
                end
                if isfield(data, 'unwrap_time_win')
                    app.unwrap_time_winEditField.Value = data.unwrap_time_win;
                else
                    app.unwrap_time_winEditField.Value = data.time_span;
                end
                if isfield(data, 'scn_time_win')
                    app.scn_time_winEditField.Value = data.scn_time_win;
                else
                    app.scn_time_winEditField.Value = data.time_span;
                end
"""
    xml = replace_once(xml, load_ui_anchor, load_ui_new, "Load callback temporal controls")

    # Remove the now-wrong coupling claim from the dataset-span label.
    xml = replace_once(
        xml,
        "app.Label_8.Text = 'Time span (note that this value will also be used for the temporal smoothing filter in StaMPS)';",
        "app.Label_8.Text = 'Dataset time span (automatically calculated from first to last acquisition)';",
        "dataset time-span label",
    )

    # These controls are source-only additions. Make the App Designer
    # round-trip limitation visible next to the existing TS Points warning.
    xml = replace_once(
        xml,
        warning_anchor,
        app_designer_warning + warning_anchor,
        "App Designer warning",
    )

    # GUI controls live on their corresponding processing-step header rows.
    weed_header = """\
            app.Step4PSweedingLabel.Position = [52 420 122 22];
            app.Step4PSweedingLabel.Text = 'Step 4 -  PS weeding';

            % Create weed_standard_devEditFieldLabel
"""
    weed_controls = """\
            app.Step4PSweedingLabel.Position = [52 420 122 22];
            app.Step4PSweedingLabel.Text = 'Step 4 -  PS weeding';

            % Create weed_time_winEditFieldLabel
            app.weed_time_winEditFieldLabel = uilabel(app.StaMPS3Tab);
            app.weed_time_winEditFieldLabel.HorizontalAlignment = 'right';
            app.weed_time_winEditFieldLabel.Position = [225 420 95 22];
            app.weed_time_winEditFieldLabel.Text = 'weed_time_win';

            % Create weed_time_winEditField
            app.weed_time_winEditField = uieditfield(app.StaMPS3Tab, 'numeric');
            app.weed_time_winEditField.Limits = [0 Inf];
            app.weed_time_winEditField.RoundFractionalValues = 'on';
            app.weed_time_winEditField.ValueDisplayFormat = '%.0f';
            app.weed_time_winEditField.FontName = 'Manrope';
            app.weed_time_winEditField.Position = [335 420 100 22];
            app.weed_time_winEditField.Value = 730;

            % Create weed_time_winDescriptionLabel
            app.weed_time_winDescriptionLabel = uilabel(app.StaMPS3Tab);
            app.weed_time_winDescriptionLabel.FontName = 'Manrope';
            app.weed_time_winDescriptionLabel.Position = [450 420 285 22];
            app.weed_time_winDescriptionLabel.Text = '(temporal smoothing window, days)';

            % Create weed_standard_devEditFieldLabel
"""
    xml = replace_once(xml, weed_header, weed_controls, "weed-time GUI header")

    unwrap_header = """\
            app.Step6PhaseunwrappingLabel.Position = [52 555 162 22];
            app.Step6PhaseunwrappingLabel.Text = 'Step 6 -  Phase unwrapping';

            % Create unwrap_grid_sizeEditFieldLabel
"""
    unwrap_controls = """\
            app.Step6PhaseunwrappingLabel.Position = [52 555 162 22];
            app.Step6PhaseunwrappingLabel.Text = 'Step 6 -  Phase unwrapping';

            % Create unwrap_time_winEditFieldLabel
            app.unwrap_time_winEditFieldLabel = uilabel(app.StaMPS4Tab);
            app.unwrap_time_winEditFieldLabel.HorizontalAlignment = 'right';
            app.unwrap_time_winEditFieldLabel.Position = [225 555 105 22];
            app.unwrap_time_winEditFieldLabel.Text = 'unwrap_time_win';

            % Create unwrap_time_winEditField
            app.unwrap_time_winEditField = uieditfield(app.StaMPS4Tab, 'numeric');
            app.unwrap_time_winEditField.Limits = [0 Inf];
            app.unwrap_time_winEditField.RoundFractionalValues = 'on';
            app.unwrap_time_winEditField.ValueDisplayFormat = '%.0f';
            app.unwrap_time_winEditField.FontName = 'Manrope';
            app.unwrap_time_winEditField.Position = [345 555 100 22];
            app.unwrap_time_winEditField.Value = 730;

            % Create unwrap_time_winDescriptionLabel
            app.unwrap_time_winDescriptionLabel = uilabel(app.StaMPS4Tab);
            app.unwrap_time_winDescriptionLabel.FontName = 'Manrope';
            app.unwrap_time_winDescriptionLabel.Position = [460 555 355 22];
            app.unwrap_time_winDescriptionLabel.Text = '(phase-noise smoothing window, days)';

            % Create unwrap_grid_sizeEditFieldLabel
"""
    xml = replace_once(xml, unwrap_header, unwrap_controls, "unwrap-time GUI header")

    scn_header = """\
            app.Step8AtmosphericfilteringLabel.Position = [50 367 178 22];
            app.Step8AtmosphericfilteringLabel.Text = 'Step 8 -  Atmospheric filtering';

            % Create scn_kriging_flagEditFieldLabel
"""
    scn_controls = """\
            app.Step8AtmosphericfilteringLabel.Position = [50 367 178 22];
            app.Step8AtmosphericfilteringLabel.Text = 'Step 8 -  Atmospheric filtering';

            % Create scn_time_winEditFieldLabel
            app.scn_time_winEditFieldLabel = uilabel(app.StaMPS5Tab);
            app.scn_time_winEditFieldLabel.HorizontalAlignment = 'right';
            app.scn_time_winEditFieldLabel.Position = [250 367 85 22];
            app.scn_time_winEditFieldLabel.Text = 'scn_time_win';

            % Create scn_time_winEditField
            app.scn_time_winEditField = uieditfield(app.StaMPS5Tab, 'numeric');
            app.scn_time_winEditField.Limits = [0 Inf];
            app.scn_time_winEditField.RoundFractionalValues = 'on';
            app.scn_time_winEditField.ValueDisplayFormat = '%.0f';
            app.scn_time_winEditField.FontName = 'Manrope';
            app.scn_time_winEditField.Position = [350 367 100 22];
            app.scn_time_winEditField.Value = 365;

            % Create scn_time_winDescriptionLabel
            app.scn_time_winDescriptionLabel = uilabel(app.StaMPS5Tab);
            app.scn_time_winDescriptionLabel.FontName = 'Manrope';
            app.scn_time_winDescriptionLabel.Position = [465 367 255 22];
            app.scn_time_winDescriptionLabel.Text = '(temporal filter window, days)';

            % Create scn_kriging_flagEditFieldLabel
"""
    xml = replace_once(xml, scn_header, scn_controls, "SCN-time GUI header")

    return xml


if __name__ == "__main__":
    if not MLAPP.is_file():
        raise SystemExit(f"Expected mlapp at {MLAPP}")
    edit_mlapp(MLAPP, patch)
    print(f"Patched {MLAPP.relative_to(REPO_ROOT)}")
