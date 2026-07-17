import zipfile


def _read_xml(path):
    with zipfile.ZipFile(path) as mlapp:
        return mlapp.read("matlab/document.xml").decode("utf-8")


def test_gacos_output_is_validated_before_stamps_step_7(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")

    generation = xml.index("aps_weather_model('gacos', 3, 3);")
    guard = xml.index("% begin TRAIN/TCA handoff guard", generation)
    step_7 = xml.index("stamps(7,7);", generation)

    assert generation < guard < step_7
    guarded = xml[guard:step_7]
    assert "fullfile(cd_fullpath, 'psver.mat')" in guarded
    assert "sprintf('tca%d.mat'" in guarded
    assert "sprintf('tca_sb%d.mat'" in guarded
    assert "who('-file', tca_path)" in guarded
    assert "PHASE_StaMPS:gacosTcaMissing" in guarded
    assert "PHASE_StaMPS:gacosTcaInvalid" in guarded


def test_resuming_step_7_rebuilds_missing_gacos_tca(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")

    resume = xml.index("% begin resumed TRAIN/TCA handoff guard")
    resumed_stamps = xml.index("stamps(stamps_first_step, stamps_last_step);", resume)
    guarded = xml[resume:resumed_stamps]

    assert "dir(fullfile(gacosFolder, '**', '*.ztd'))" in guarded
    assert "setparm_aps('gacos_datapath', gacosFolder)" in guarded
    assert "aps_weather_model('gacos', 3, 3);" in guarded
    assert "who('-file', tca_path)" in guarded
    assert "PHASE_StaMPS:gacosMapsMissing" in guarded


def test_gacos_generation_is_pinned_to_processing_root(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    call = "aps_weather_model('gacos', 3, 3);"

    for offset in (i for i in range(len(xml)) if xml.startswith(call, i)):
        context = xml[max(0, offset - 180):offset + len(call) + 180]
        assert "cd(cd_fullpath);" in context
