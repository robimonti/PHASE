import zipfile


def _read_xml(path):
    with zipfile.ZipFile(path) as mlapp:
        return mlapp.read("matlab/document.xml").decode("utf-8")


def _start_callback(xml):
    start = xml.index("function StartButtonPushed(app, event)")
    end = xml.index("% Button pushed function: LoadButton", start)
    return xml[start:end]


def test_train_request_requires_both_user_controls(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)

    assert (
        "train_requested = (train_flag == 0) && "
        "strcmpi(strtrim(subtr_tropo), 'y');"
    ) in callback


def test_one_effective_flag_controls_processing_and_export(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)

    decision = callback.index("% begin effective tropospheric-correction decision")
    processing = callback.index("%% -------------------------- RUN STAMPS")
    assert decision < processing
    assert (
        "tropo_correction_enabled = (train_flag == 0) && ..."
    ) in callback[decision:processing]

    # Two preparation branches, StaMPS processing and displacement export use
    # the flag directly; atmosphere export combines it with other predicates.
    assert callback.count("if tropo_correction_enabled") == 4
    assert (
        "stamps_last_step > 6 && tropo_correction_enabled && "
        "contains(ph_output, 'unwrapped')"
    ) in callback
    assert "if train_flag == 0 % TRAIN tropo correction included" not in callback


def test_subtr_tropo_n_run_does_not_enter_tca2_export(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)
    atmosphere = callback.index("%% ------------------ EXPORT ATMOSPHERIC CORRECTION")
    tca_load = callback.index("load('tca2.mat'", atmosphere)
    guard = callback[atmosphere:tca_load]

    assert "tropo_correction_enabled" in guard
    assert "train_flag == 0" not in guard
