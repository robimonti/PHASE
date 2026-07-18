"""Stop PHASE_StaMPS when visible tropo controls were not saved.

Start historically loaded input_StaMPS.mat without warning when controls in the
GUI had been edited.  Silently auto-saving every control is unsafe because the
app does not auto-load the MAT file at startup: untouched controls may still
contain compiled defaults.  Compare only the three controls that decide TRAIN
and GACOS behavior, and require the explicit Load/Edit/Save workflow on a
mismatch without overwriting the existing configuration.
"""

from pathlib import Path

from mlapp_roundtrip import edit_mlapp


REPO_ROOT = Path(__file__).resolve().parents[1]
MLAPP = REPO_ROOT / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"


UNSAFE_SNAPSHOT = """\
                % begin Start/current-UI configuration snapshot
                % Start must execute the values currently visible in the app.
                % Previously it loaded input_StaMPS.mat directly, so an edited
                % control (notably subtr_tropo='n') silently kept its old saved
                % value unless the separate Save button had been pressed first.
                SaveButtonPushed(app, event);
                if ~isequal(app.StatusLamp.Color, [0, 1, 0])
                    error('PHASE_StaMPS:configurationSaveFailed', ...
                        ['Could not save the visible settings to input_StaMPS.mat. ', ...
                         'Processing was not started.']);
                end
                updateOutput(app, ['Using visible configuration: subtr_tropo=' ...
                    char(string(app.subtr_tropoEditField.Value)) ...
                    ', tropo_method=' char(string(app.tropo_methodEditField.Value))]);
                % end Start/current-UI configuration snapshot

"""


LOAD_END = """\
                    'scla_method', 'scla_drop_index', 'scn_wavelength', 'scn_kriging_flag', 'ref_centre_lonlat', ...
                    'ref_radius', 'ref_velocity', 'plot_s', 'ref_centre_lonlat_w', 'ref_radius_w', 'ph_output');

                % AutoDetectParameters (startupFcn) reads ground-truth values
"""


LOAD_END_WITH_GUARD = """\
                    'scla_method', 'scla_drop_index', 'scn_wavelength', 'scn_kriging_flag', 'ref_centre_lonlat', ...
                    'ref_radius', 'ref_velocity', 'plot_s', 'ref_centre_lonlat_w', 'ref_radius_w', 'ph_output');

                % begin saved/visible tropospheric-configuration guard
                % The app deliberately has separate Load and Save buttons. If
                % a user edits these controls without Save, Start would load an
                % older value and could unexpectedly launch GACOS or seek tca2.
                % Never overwrite the MAT file here: other untouched controls
                % may still contain the app's compiled defaults until Load.
                loaded_subtr_tropo = char(string(subtr_tropo));
                visible_subtr_tropo = char(string(app.subtr_tropoEditField.Value));
                loaded_tropo_method = char(string(tropo_method));
                visible_tropo_method = char(string(app.tropo_methodEditField.Value));
                visible_train_flag = double(~app.TRAINatmosphericcorrectionCheckBox.Value);
                tropo_config_mismatch = ...
                    ~strcmpi(strtrim(loaded_subtr_tropo), strtrim(visible_subtr_tropo)) || ...
                    ~strcmpi(strtrim(loaded_tropo_method), strtrim(visible_tropo_method)) || ...
                    train_flag ~= visible_train_flag;
                if tropo_config_mismatch
                    error('PHASE_StaMPS:unsavedConfiguration', ...
                        ['The visible TRAIN/tropospheric controls do not match input_StaMPS.mat. ', ...
                         'Processing was not started and the MAT file was not changed. ', ...
                         'Press Load, make the intended changes, press Save, then press Start. ', ...
                         'Saved: TRAIN=%d, subtr_tropo=%s, tropo_method=%s. ', ...
                         'Visible: TRAIN=%d, subtr_tropo=%s, tropo_method=%s.'], ...
                         train_flag == 0, loaded_subtr_tropo, loaded_tropo_method, ...
                         visible_train_flag == 0, visible_subtr_tropo, visible_tropo_method);
                end
                updateOutput(app, sprintf( ...
                    'Loaded configuration: TRAIN=%d, subtr_tropo=%s, tropo_method=%s', ...
                    train_flag == 0, loaded_subtr_tropo, loaded_tropo_method));
                clear loaded_subtr_tropo visible_subtr_tropo loaded_tropo_method
                clear visible_tropo_method visible_train_flag tropo_config_mismatch
                % end saved/visible tropospheric-configuration guard

                % AutoDetectParameters (startupFcn) reads ground-truth values
"""


def patch(xml: str) -> str:
    marker = "% begin saved/visible tropospheric-configuration guard"
    if marker in xml:
        print("PHASE_StaMPS.mlapp already guards unsaved tropo controls.")
        return xml

    if UNSAFE_SNAPSHOT in xml:
        xml = xml.replace(UNSAFE_SNAPSHOT, "", 1)
    elif "% begin Start/current-UI configuration snapshot" in xml:
        raise RuntimeError("Unrecognized Start auto-save draft; refusing broad replacement")

    if xml.count(LOAD_END) != 1:
        raise RuntimeError("input_StaMPS.mat load-end anchor was not found exactly once")
    return xml.replace(LOAD_END, LOAD_END_WITH_GUARD, 1)


if __name__ == "__main__":
    if not MLAPP.is_file():
        raise SystemExit(f"Expected mlapp at {MLAPP}")
    edit_mlapp(MLAPP, patch)
    print(f"Patched {MLAPP.relative_to(REPO_ROOT)}")
