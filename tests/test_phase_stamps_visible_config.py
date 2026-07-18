import zipfile


def _read_xml(path):
    with zipfile.ZipFile(path) as mlapp:
        return mlapp.read("matlab/document.xml").decode("utf-8")


def _start_callback(xml):
    start = xml.index("function StartButtonPushed(app, event)")
    end = xml.index("% Button pushed function: LoadButton", start)
    return xml[start:end]


def test_start_checks_visible_tropo_controls_after_loading_mat(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)

    load_mat = callback.index("load('input_StaMPS.mat'")
    guard = callback.index("% begin saved/visible tropospheric-configuration guard")
    train_decision = callback.index("train_requested =", guard)
    assert load_mat < guard < train_decision

    guarded = callback[guard:train_decision]
    assert "app.subtr_tropoEditField.Value" in guarded
    assert "app.tropo_methodEditField.Value" in guarded
    assert "app.TRAINatmosphericcorrectionCheckBox.Value" in guarded
    assert "PHASE_StaMPS:unsavedConfiguration" in guarded


def test_start_never_overwrites_mat_to_resolve_a_visible_mismatch(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)

    assert "SaveButtonPushed(app, event);" not in callback
    assert "Processing was not started and the MAT file was not changed" in callback
    assert "Press Load, make the intended changes, press Save, then press Start" in callback


def test_start_reports_the_loaded_effective_tropo_controls(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)

    assert "Loaded configuration: TRAIN=%d, subtr_tropo=%s, tropo_method=%s" in callback
