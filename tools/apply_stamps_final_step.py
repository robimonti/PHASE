"""Make PHASE_StaMPS honor the selected final StaMPS step.

The no-atmospheric-correction branch historically ran to a literal Step 7
and then repeated processing from Step 6.  Replace that sequence with one
call using the saved first/last-step range.  Also include both selectors in
the saved-versus-visible guard so an unsaved GUI value cannot be mistaken for
the value actually loaded from input_StaMPS.mat.
"""

from pathlib import Path

from mlapp_roundtrip import edit_mlapp


REPO_ROOT = Path(__file__).resolve().parents[1]
MLAPP = REPO_ROOT / "legacy" / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"


def replace_once(xml: str, old: str, new: str, description: str) -> str:
    if xml.count(old) != 1:
        raise RuntimeError(f"{description} anchor was not found exactly once")
    return xml.replace(old, new, 1)


def patch(xml: str) -> str:
    marker = "% begin selected StaMPS step guard"
    if marker in xml:
        print("PHASE_StaMPS.mlapp already honors the selected final step.")
        return xml

    visible_anchor = """\
                visible_scn_time_win = app.scn_time_winEditField.Value;
                tropo_config_mismatch = ...
"""
    visible_new = """\
                visible_scn_time_win = app.scn_time_winEditField.Value;
                % begin selected StaMPS step guard
                loaded_stamps_first_step = str2double(char(string(stamps_first_step)));
                loaded_stamps_last_step = str2double(char(string(stamps_last_step)));
                visible_stamps_first_step = str2double(char(string(app.StaMPSfirststepDropDown.Value)));
                visible_stamps_last_step = str2double(char(string(app.StaMPSlaststepDropDown.Value)));
                % end selected StaMPS step guard
                tropo_config_mismatch = ...
"""
    xml = replace_once(xml, visible_anchor, visible_new, "visible step values")

    mismatch_anchor = """\
                time_window_config_mismatch = ...
                    ~isequaln(weed_time_win, visible_weed_time_win) || ...
                    ~isequaln(unwrap_time_win, visible_unwrap_time_win) || ...
                    ~isequaln(scn_time_win, visible_scn_time_win);
                if tropo_config_mismatch || time_window_config_mismatch
"""
    mismatch_new = """\
                time_window_config_mismatch = ...
                    ~isequaln(weed_time_win, visible_weed_time_win) || ...
                    ~isequaln(unwrap_time_win, visible_unwrap_time_win) || ...
                    ~isequaln(scn_time_win, visible_scn_time_win);
                step_range_config_mismatch = ...
                    ~isequaln(loaded_stamps_first_step, visible_stamps_first_step) || ...
                    ~isequaln(loaded_stamps_last_step, visible_stamps_last_step);
                if tropo_config_mismatch || time_window_config_mismatch || step_range_config_mismatch
"""
    xml = replace_once(xml, mismatch_anchor, mismatch_new, "step mismatch condition")

    error_anchor = """\
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
    error_new = """\
                        ['The visible TRAIN/tropospheric/temporal-window/step controls do not match input_StaMPS.mat. ', ...
                         'Processing was not started and the MAT file was not changed. ', ...
                         'Press Load, make the intended changes, press Save, then press Start. ', ...
                         'Saved: TRAIN=%d, subtr_tropo=%s, tropo_method=%s, time windows=%g/%g/%g days, steps=%g->%g. ', ...
                         'Visible: TRAIN=%d, subtr_tropo=%s, tropo_method=%s, time windows=%g/%g/%g days, steps=%g->%g.'], ...
                         train_flag == 0, loaded_subtr_tropo, loaded_tropo_method, ...
                         weed_time_win, unwrap_time_win, scn_time_win, ...
                         loaded_stamps_first_step, loaded_stamps_last_step, ...
                         visible_train_flag == 0, visible_subtr_tropo, visible_tropo_method, ...
                         visible_weed_time_win, visible_unwrap_time_win, visible_scn_time_win, ...
                         visible_stamps_first_step, visible_stamps_last_step);
"""
    xml = replace_once(xml, error_anchor, error_new, "unsaved step error")

    status_anchor = """\
                    ['Loaded configuration: TRAIN=%d, subtr_tropo=%s, tropo_method=%s, ', ...
                     'weed/unwrap/scn time windows=%g/%g/%g days'], ...
                    train_flag == 0, loaded_subtr_tropo, loaded_tropo_method, ...
                    weed_time_win, unwrap_time_win, scn_time_win));
"""
    status_new = """\
                    ['Loaded configuration: TRAIN=%d, subtr_tropo=%s, tropo_method=%s, ', ...
                     'weed/unwrap/scn time windows=%g/%g/%g days, StaMPS steps=%g->%g'], ...
                    train_flag == 0, loaded_subtr_tropo, loaded_tropo_method, ...
                    weed_time_win, unwrap_time_win, scn_time_win, ...
                    loaded_stamps_first_step, loaded_stamps_last_step));
"""
    xml = replace_once(xml, status_anchor, status_new, "loaded step status")

    clear_anchor = """\
                clear time_window_config_mismatch legacy_time_window_fields
                % end saved/visible tropospheric-configuration guard
"""
    clear_new = """\
                clear time_window_config_mismatch legacy_time_window_fields
                clear loaded_stamps_first_step loaded_stamps_last_step
                clear visible_stamps_first_step visible_stamps_last_step
                clear step_range_config_mismatch
                % end saved/visible tropospheric-configuration guard
"""
    xml = replace_once(xml, clear_anchor, clear_new, "step guard cleanup")

    no_tropo_anchor = """\
                        setparm('subtr_tropo', 'n');
                        stamps(stamps_first_step,7); % execute StaMPS steps from first step to 7
                        stamps(6,stamps_last_step); % subtract the computed corrections before the phase unwrapping
"""
    no_tropo_new = """\
                        setparm('subtr_tropo', 'n');
                        stamps(stamps_first_step, stamps_last_step); % execute the complete range selected in the GUI
"""
    xml = replace_once(xml, no_tropo_anchor, no_tropo_new, "no-tropo StaMPS range")

    return xml


if __name__ == "__main__":
    if not MLAPP.is_file():
        raise SystemExit(f"Expected mlapp at {MLAPP}")
    edit_mlapp(MLAPP, patch)
    print(f"Patched {MLAPP.relative_to(REPO_ROOT)}")
