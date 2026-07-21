import zipfile


def _read_xml(path):
    with zipfile.ZipFile(path) as mlapp:
        return mlapp.read("matlab/document.xml").decode("utf-8")


def test_temporal_windows_have_independent_gui_controls_and_defaults(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")

    expected = {
        "weed_time_win": ("StaMPS3Tab", "730"),
        "unwrap_time_win": ("StaMPS4Tab", "730"),
        "scn_time_win": ("StaMPS5Tab", "365"),
    }
    for name, (tab, default) in expected.items():
        assert f"{name}EditField" in xml
        declaration = xml.index(f"{name}EditField")
        assert "matlab.ui.control.NumericEditField" in xml[declaration:declaration + 100]
        assert f"app.{name}EditField = uieditfield(app.{tab}, 'numeric');" in xml
        assert f"app.{name}EditField.Value = {default};" in xml


def test_save_start_and_load_persist_all_three_windows(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")

    for name in ("weed_time_win", "unwrap_time_win", "scn_time_win"):
        assert f"{name} = app.{name}EditField.Value;" in xml
        assert f"'{name}'" in xml
        assert f"isfield(data, '{name}')" in xml
        assert f"app.{name}EditField.Value = data.{name};" in xml

    assert xml.count(
        "'time_span', 'weed_time_win', 'unwrap_time_win', 'scn_time_win', ..."
    ) == 2


def test_old_mat_files_fall_back_to_dataset_span(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    start = xml.index("% begin independent StaMPS temporal windows")
    end = xml.index("% end independent StaMPS temporal windows", start)
    fallback = xml[start:end]

    for name in ("weed_time_win", "unwrap_time_win", "scn_time_win"):
        assert f"exist('{name}', 'var') ~= 1" in fallback
        assert f"{name} = time_span;" in fallback
        assert f"app.{name}EditField.Value = {name};" not in fallback
    assert "Press Save to migrate the MAT file" in fallback


def test_runtime_no_longer_uses_dataset_span_as_processing_window(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")

    for name in ("weed_time_win", "unwrap_time_win", "scn_time_win"):
        assert f"setparm('{name}', {name});" in xml
        assert f"setparm('{name}', time_span);" not in xml

    assert "Dataset time span (automatically calculated" in xml
    assert "Time span (note that this value will also be used" not in xml


def test_unsaved_time_window_edits_are_guarded(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    guard_start = xml.index("% begin saved/visible tropospheric-configuration guard")
    guard_end = xml.index("% end saved/visible tropospheric-configuration guard")
    guard = xml[guard_start:guard_end]

    assert "time_window_config_mismatch" in guard
    assert "tropo_config_mismatch || time_window_config_mismatch" in guard
    for name in ("weed_time_win", "unwrap_time_win", "scn_time_win"):
        assert f"app.{name}EditField.Value" in guard


def test_installer_seeds_independent_time_window_defaults(phase_root):
    installer = (phase_root / "installer" / "install-phase.ps1").read_text(
        encoding="utf-8"
    )

    assert '"    weed_time_win = 730;"' in installer
    assert '"    unwrap_time_win = 730;"' in installer
    assert '"    scn_time_win = 365;"' in installer
    assert (
        "'time_span', 'weed_time_win', 'unwrap_time_win', 'scn_time_win', "
        "'year_0'"
    ) in installer
