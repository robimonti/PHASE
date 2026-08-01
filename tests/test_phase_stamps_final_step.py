import zipfile


def _read_xml(path):
    with zipfile.ZipFile(path) as mlapp:
        return mlapp.read("matlab/document.xml").decode("utf-8")


def _start_callback(xml):
    start = xml.index("function StartButtonPushed(app, event)")
    end = xml.index("% Button pushed function: LoadButton", start)
    return xml[start:end]


def _no_tropo_branch(callback):
    start = callback.index("% begin TRAIN-degradation GUI notice")
    end = callback.index("elseif contains(ph_output, 'wrapped')", start)
    return callback[start:end]


def test_no_tropo_processing_uses_selected_full_range_once(phase_root):
    xml = _read_xml(phase_root / "legacy" / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    branch = _no_tropo_branch(_start_callback(xml))

    assert "stamps(stamps_first_step, stamps_last_step);" in branch
    assert "stamps(stamps_first_step,7)" not in branch
    assert "stamps(6,stamps_last_step)" not in branch
    assert branch.count("stamps(") == 1


def test_start_rejects_unsaved_first_or_last_step(phase_root):
    xml = _read_xml(phase_root / "legacy" / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)
    guard_start = callback.index("% begin saved/visible tropospheric-configuration guard")
    guard_end = callback.index("% end saved/visible tropospheric-configuration guard")
    guard = callback[guard_start:guard_end]

    assert "app.StaMPSfirststepDropDown.Value" in guard
    assert "app.StaMPSlaststepDropDown.Value" in guard
    assert "loaded_stamps_last_step" in guard
    assert "visible_stamps_last_step" in guard
    assert "step_range_config_mismatch" in guard
    assert "steps=%g->%g" in guard


def test_loaded_status_reports_effective_step_range(phase_root):
    xml = _read_xml(phase_root / "legacy" / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)

    assert "StaMPS steps=%g->%g" in callback


def test_train_gacos_intermediate_step_7_is_preserved(phase_root):
    xml = _read_xml(phase_root / "legacy" / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    callback = _start_callback(xml)

    # This is not a final-step limit: TRAIN needs the Step-7 product before
    # rerunning Step 6 through the user-selected last step.
    assert "stamps(7,7); % execute StaMPS step 7" in callback
    assert "stamps(6,stamps_last_step); % subtract the computed corrections" in callback
