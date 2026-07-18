def test_diagnostic_captures_step6_runtime_evidence(phase_root):
    source = (
        phase_root / "PHASE_Preprocessing" / "diagnose_PHASE_StaMPS.m"
    ).read_text(encoding="utf-8")

    assert "input_cfg = load(input_file, 'installation_folder', 'train_flag'" in source
    assert "Saved subtr_tropo:" in source
    assert "Saved tropo_method:" in source
    assert "Saved TRAIN enabled:" in source
    assert "external', 'snaphu', 'bin', 'snaphu.exe" in source
    assert "which('ps_unwrap', '-all')" in source
    assert "which('uw_stat_costs', '-all')" in source
    assert "uw_nosnaphu" in source
    assert "system('where snaphu')" in source
    assert "PHASE_StaMPS_diagnostic.txt" in source
