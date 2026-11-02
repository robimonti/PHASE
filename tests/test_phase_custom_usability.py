"""Interface contract for the streamlined custom modelling workflow."""


def test_custom_controls_use_ram_and_no_separate_time_offset(phase_root):
    schema = (phase_root / "+phase_model_beta/schema.m").read_text()
    defaults = (phase_root / "+phase_model_beta/defaultConfig.m").read_text()
    assert "{'Automatic','Custom'}" in schema
    assert "customTimeOffset" not in schema + defaults
    assert "customMaxObservations" not in schema + defaults
    assert "customMemoryGB" in schema and "'GiB'" in schema
    assert "'days')" in schema.split("field('step_t_ST'", 1)[1].split("field(", 1)[0]


def test_custom_prediction_reports_to_host_and_avoids_dense_repeated_products(phase_root):
    run = (phase_root / "+phase_model_custom/run.m").read_text()
    predict = (phase_root / "+phase_model_custom/predict.m").read_text()
    engine = (phase_root / "+phase_model_beta/ProcessingEngine.m").read_text()
    assert "customConfig.progressCallback" in engine
    assert "Custom prediction:" not in run
    for phase in ("Output grid", "Input PS", "Query points"):
        assert phase in run
    assert "V=residual.factor\\K'" in predict
    assert "base=targetVariance-sum(V.^2,1)'" in predict


def test_private_student_input_excluded(phase_root):
    import importlib.util
    import subprocess

    result = subprocess.run(
        ["git", "check-ignore", "examples/selected_random_points_DSC_VV.xlsx"],
        cwd=phase_root, capture_output=True, text=True,
    )
    assert result.returncode == 0
    spec = importlib.util.spec_from_file_location(
        "phase_private_input_check", phase_root / "installer/install-phase-unix.py"
    )
    installer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(installer)
    installer.SOURCE_FOR_COPY = phase_root
    assert installer.ignored(str(phase_root / "examples"), [
        "selected_random_points_DSC_VV.xlsx", "custom_models"
    ]) == {"selected_random_points_DSC_VV.xlsx"}


def test_coordinate_overrides_are_optional_and_saved_origins_are_expanded(phase_root):
    ui = (phase_root / "phase_model_beta_ui/app.js").read_text()
    assert "Coordinate settings · optional" in ui
    assert 'details.open = c.customCoordinates === "lonlat"' in ui
    assert "explicitOrigin(c.customOriginX)" in ui
    assert "explicitOrigin(c.customOriginY)" in ui
