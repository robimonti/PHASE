def test_manual_points_materialize_through_existing_file_path(phase_root):
    defaults=(phase_root / "+phase_model_beta/defaultConfig.m").read_text()
    schema=(phase_root / "+phase_model_beta/schema.m").read_text()
    assert "'queryPointSource','file'" in defaults
    assert "'queryPoints','[]'" in defaults
    assert "'queryPointSource'" in schema and "'queryPoints'" in schema
    save=(phase_root / "+phase_model_beta/saveConfig.m").read_text()
    assert "materializeQueryPoints(config,rootDir)" in save
    writer=(phase_root / "+phase_model_beta/materializeQueryPoints.m").read_text()
    assert "query_points_gui.txt" in writer and "parseQueryPoints" in writer
    ui=(phase_root / "phase_model_beta_ui/app.js").read_text()
    assert "createQueryPointEditor" in ui and "+ Add point" in ui


def test_covariance_diagnostics_export_actual_kernel_and_counts(phase_root):
    for adapter in ("nativeCovariance.m","nativeCovariance1D.m"):
        code=(phase_root / "+phase_model_custom" / adapter).read_text()
        assert "settings.pairCountTime" in code and "settings.pairCountSpace" in code
        assert "settings.marginalTimeParameters" in code
    prepare=(phase_root / "+phase_model_custom/prepareResidual.m").read_text()
    assert "exportCovarianceDiagnostics" in prepare
    export=(phase_root / "+phase_model_custom/exportCovarianceDiagnostics.m").read_text()
    assert "savefig(" in export and "print(f," in export
    assert "FittedSignal_mm2" in export and "PairCount" in export
    assert "subplot(1,2,1)" in export and "subplot(1,2,2)" in export
    assert "legend('empirical', 'signal', 'noise')" in export
