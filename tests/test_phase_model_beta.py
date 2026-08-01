import importlib.util
import re
import zipfile


def _text(path):
    return path.read_text(encoding="utf-8")


def _stable_xml(phase_root):
    with zipfile.ZipFile(phase_root / "legacy/PHASE_model.mlapp") as app:
        return app.read("matlab/document.xml").decode("utf-8")


def test_model_beta_has_standalone_text_launcher(phase_root):
    launcher = _text(phase_root / "PHASE_Model_beta.m")
    assert "phase_model_beta.App" in launcher
    assert "PHASE_model.mlapp" not in launcher
    assert not (phase_root / "PHASE_Model_beta.mlapp").exists()


def test_model_beta_engine_is_reproducibly_extracted(phase_root):
    tool_path = phase_root / "tools" / "extract_model_beta_engine.py"
    spec = importlib.util.spec_from_file_location("extract_model_beta_engine", tool_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    generated = _text(phase_root / "+phase_model_beta" / "LegacyEngine.m")
    assert generated == module.extract(_stable_xml(phase_root))


def test_model_beta_preserves_every_stable_callback(phase_root):
    stable = _stable_xml(phase_root)
    engine = _text(phase_root / "+phase_model_beta" / "LegacyEngine.m")
    stable_functions = re.findall(
        r"^\s*function\s+(?:\[[^]]+\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(",
        stable,
        re.M,
    )
    engine_functions = re.findall(
        r"^\s*function\s+(?:\[[^]]+\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(",
        engine,
        re.M,
    )
    assert len(engine_functions) == len(stable_functions) + 1
    assert set(stable_functions).issubset(engine_functions)
    assert "notifyBetaProgress" in engine_functions
    assert "StartButtonPushed" in engine_functions
    assert "SaveButtonPushed" in engine_functions
    assert "LoadButtonPushed" in engine_functions
    assert "phase_model_beta.projectRoot()" in engine
    assert "phase_model_beta.themeLegacyEngine(app)" in engine
    assert "phase_model_beta.resolvePythonPath(app.pythonPath)" in engine
    assert "app.UIFigure.Visible = 'off'" in engine
    assert "phase_model_beta.runCommandHidden(app, command" in engine
    assert "waitbar(" not in engine
    assert "mfilename('fullpath')" not in engine


def test_model_beta_reuses_installer_python_and_has_self_test(phase_root):
    resolver = _text(phase_root / "+phase_model_beta" / "resolvePythonPath.m")
    self_test = _text(phase_root / "+phase_model_beta" / "selfTest.m")

    assert "APPDATA" in resolver
    assert "'PHASE', 'python.txt'" in resolver
    assert "input_model.mat remains authoritative" in resolver
    assert "PHASE Model beta standalone self-test passed." in self_test


def test_model_beta_has_modern_html_shell_and_live_monitor(phase_root):
    controller = _text(phase_root / "+phase_model_beta" / "App.m")
    html = _text(phase_root / "phase_model_beta_ui" / "index.html")
    js = _text(phase_root / "phase_model_beta_ui" / "app.js")
    css = _text(phase_root / "phase_model_beta_ui" / "styles.css")

    assert "phase_model_beta.LegacyEngine()" in controller
    assert "obj.Engine.UIFigure.Visible = 'off'" in controller
    assert "diary on" in controller
    assert "liveLogUrl" in controller
    assert "Run monitor" in html
    assert "Module 2" in html
    assert "Module 2 · Beta" not in html
    assert "type = \"date\"" in js
    assert "pollLiveLog" in js
    assert 'id="stop-button"' in html
    assert 'send("Stop", {})' in js
    assert "StopRequested" in controller
    assert "--blue: rgb(53, 101, 207)" in css


def test_model_beta_save_uses_correct_engine_callback_arity(phase_root):
    controller = _text(phase_root / "+phase_model_beta" / "App.m")

    assert "LoadButtonPushed([],[])" not in controller
    assert "StartButtonPushed([],[])" not in controller
    assert "LoadButtonPushed([])" not in controller
    assert controller.count("phase_model_beta.applyConfigToEngine") == 4
    assert "StartButtonPushed([])" in controller


def test_model_beta_dynamic_sections_and_pair_layout(phase_root):
    schema = _text(phase_root / "+phase_model_beta" / "schema.m")
    js = _text(phase_root / "phase_model_beta_ui" / "app.js")
    css = _text(phase_root / "phase_model_beta_ui" / "styles.css")

    assert 'group.id === "observations" || group.id === "splines"' in js
    assert 'group.id === "spatial1d"' in js
    assert 'group.id === "spatial2d"' in js
    assert "pairRole(meta.id)" in js
    assert ".field-card.pair-method" in css
    assert ".field-card.pair-value" in css
    assert "covarianceModelRole(meta.id)" in js
    assert ".field-card.covariance-model" in css
    assert "spline_method_DET2D: deterministic &&" in js
    assert "Export interpolated PS observations" in schema


def test_model_beta_engine_does_not_load_config_through_hidden_widgets(phase_root):
    package = phase_root / "+phase_model_beta"
    controller = _text(package / "App.m")
    engine = _text(package / "LegacyEngine.m")
    applicator = _text(package / "applyConfigToEngine.m")

    assert "Configuration is applied by the standalone controller." in engine
    assert "% Automatic load from input_model.mat if it exists" not in engine
    assert "engine.(name) = config.(name)" in applicator
    assert "phase_model_beta.applyConfigToEngine(obj.Engine,candidate)" in controller

    start = engine[
        engine.index("function StartButtonPushed(app, event)") :
        engine.index("% Value changed function: lonmaxEditField")
    ]
    for hidden_widget in (
        "app.lonminEditField.Value",
        "app.MethodDropDown_temporal.Value",
        "app.manualvalueEditField_noise.Value",
        "app.xnEditField.Value",
    ):
        assert hidden_widget not in start
    assert "lonMinAOI = app.lonMinAOI" in start
    assert "dtCov_STC2D = app.dtCov_STC2D" in start


def test_model_beta_can_estimate_full_ps_extent(phase_root):
    package = phase_root / "+phase_model_beta"
    controller = _text(package / "App.m")
    estimator = _text(package / "estimatePsBoundingBox.m")
    js = _text(phase_root / "phase_model_beta_ui" / "app.js")

    assert "case 'estimatepsbounds'" in controller
    assert "phase_model_beta.estimatePsBoundingBox" in controller
    assert "values(3:end,2:3)" in estimator
    assert 'send("EstimatePsBounds", payload())' in js
    assert "Select full PS extent" in js


def test_model_beta_initialises_map_from_ps_extent(phase_root):
    controller = _text(phase_root / "+phase_model_beta" / "App.m")

    assert "obj.applyDefaultPsBounds(false);" in controller
    assert "obj.applyDefaultPsBounds(true);" in controller
    assert "Default AOI fitted to the PS extent" in controller
    assert "obj.Config.aoi_polygon_lonlat = bboxPolygon(bounds)" in controller
    assert "obj.MapPsExtent = bboxPolygon(bounds)" in controller
    assert "'id','ps-extent','name','PS extent'" in controller
    assert "footprints," in _text(
        phase_root / "phase_model_beta_ui" / "app.js"
    )


def test_model_beta_displays_and_uses_the_same_shapefile_aoi(phase_root):
    package = phase_root / "+phase_model_beta"
    controller = _text(package / "App.m")
    engine = _text(package / "LegacyEngine.m")
    reader = _text(package / "readAoiShapefile.m")
    js = _text(phase_root / "phase_model_beta_ui" / "app.js")

    assert controller.count("phase_model_beta.readAoiShapefile") >= 1
    assert "obj.MapAoiFootprints = features" in controller
    assert "phase_model_beta.readAoiShapefile" in engine
    assert "finitePair = isfinite(x) & isfinite(y)" in reader
    assert "lonlat(end+1,:) = [NaN NaN]" in reader
    assert "selectedFootprints.flatMap" in js
    assert "fitToFootprints(false)" in js


def test_model_beta_rejects_empty_aoi_and_never_reverse_geocodes_nan(phase_root):
    engine = _text(phase_root / "+phase_model_beta" / "LegacyEngine.m")
    reverse = _text(phase_root / "MatlabFunctions" / "get_place_from_coordinates.m")

    assert "PHASE_Model_beta:noPsInsideAoi" in engine
    assert "isfinite(query_lon) && isfinite(query_lat)" in engine
    assert "~isfinite(lon) || ~isfinite(lat)" in reverse
    assert "Reverse geocoding skipped" in reverse


def test_model_beta_uses_stable_output_root_and_direct_geosplinter_runner(phase_root):
    engine = _text(phase_root / "+phase_model_beta" / "LegacyEngine.m")
    runner = _text(phase_root / "MatlabFunctions" / "runGeoSplinter.m")
    model_files = [
        phase_root / "MatlabFunctions" / "ModellingInTime.m",
        phase_root / "MatlabFunctions" / "STmodel_DET1D.m",
        phase_root / "MatlabFunctions" / "STmodel_DET2D.m",
        phase_root / "MatlabFunctions" / "STmodel_STC1D.m",
        phase_root / "MatlabFunctions" / "STmodel_STC2D.m",
    ]

    assert "outputRoot = fileparts(runtimeRoot)" in engine
    assert "outputDir = fullfile('..',outputDir)" in engine
    assert "app.outputDir = char(java.io.File(outputDir).getCanonicalPath())" in engine
    assert "redirectInput(java.io.File(jobFile))" in runner
    assert "endsWith(lower(executable),'.exe')" in runner
    combined = "\n".join(_text(path) for path in model_files)
    assert combined.count("runGeoSplinter(") == 13
    assert "temp_bat" not in combined
    assert "system(job_execution" not in combined

    smoke_test = _text(
        phase_root / "+phase_model_beta" / "geoSplinterSelfTest.m"
    )
    assert "runGeoSplinter(executable,jobPath)" in smoke_test
    assert "outputName = 'PS_1_cub'" in smoke_test
    assert "'.par.txt','.std.txt','.mat.txt'" in smoke_test
    assert "PHASE Model geoSplinter native-runtime self-test passed." in smoke_test


def test_temporal_interpolation_consolidates_duplicate_sample_points(phase_root):
    helper = _text(phase_root / "MatlabFunctions" / "interp1Unique.m")
    modelling = _text(phase_root / "MatlabFunctions" / "ModellingInTime.m")
    fitters = "\n".join(
        _text(phase_root / "MatlabFunctions" / name)
        for name in (
            "fit_exp_auto.m",
            "fit_gaussian_auto.m",
            "fit_gaussian_cosine_auto.m",
        )
    )

    assert "[uniqueX,~,groups] = unique(x)" in helper
    assert "accumarray(groups,y" in helper
    assert "values = interp1(uniqueX,uniqueY" in helper
    assert modelling.count("interp1Unique(") == 3
    assert fitters.count("interp1Unique(") == 6
    assert "obj.appendLog(getReport(ME,'extended','hyperlinks','off'))" in _text(
        phase_root / "+phase_model_beta" / "App.m"
    )


def test_model_beta_supports_interactive_polygon_aoi(phase_root):
    package = phase_root / "+phase_model_beta"
    defaults = _text(package / "defaultConfig.m")
    validator = _text(package / "validateConfig.m")
    engine = _text(package / "LegacyEngine.m")
    controller = _text(package / "App.m")
    html = _text(phase_root / "phase_model_beta_ui" / "index.html")
    js = _text(phase_root / "phase_model_beta_ui" / "app.js")

    assert "aoi_polygon_lonlat" in defaults
    assert "polyarea" in validator
    assert "lonlatAOI = app.aoi_polygon_lonlat" in engine
    assert "case 'mapaoichanged'" in controller
    assert 'id="aoi-map-panel"' in html
    assert "new PhaseMap" in js
    assert 'send("MapAoiChanged"' in js
    assert '<script src="map.js"></script>' in html
    assert "../PHASE_Preprocessing/phase_preprocessing_beta_ui/map.js" not in html
    assert "copyfile(mapSource,mapTarget,'f')" in controller
    assert "obj.ensureAssets();" in controller


def test_model_beta_temporal_mode_skips_spatial_grid_and_exports_report_figures(phase_root):
    engine = _text(phase_root / "+phase_model_beta" / "LegacyEngine.m")
    exporter = _text(phase_root / "+phase_model_beta" / "exportFigure.m")

    assert "if strcmp(procType, 'temporal')" in engine
    assert "Pure temporal mode: spatial grid generation skipped." in engine
    assert engine.count("phase_model_beta.exportFigure") == 5
    assert "exportgraphics" in exporter
    assert "figureExportFailed" in exporter


def test_model_beta_hard_stop_reaches_temporal_loop_and_external_helpers(phase_root):
    package = phase_root / "+phase_model_beta"
    engine = _text(package / "LegacyEngine.m")
    modelling = _text(phase_root / "MatlabFunctions" / "ModellingInTime.m")
    runner = _text(package / "runCommandHidden.m")

    assert "phase_model_beta.throwIfStopped(app)" in engine
    assert "'stop_check'" in engine
    assert "addParameter(p, 'stop_check'" in modelling
    assert "stop_check();" in modelling
    assert "process.destroyForcibly()" in runner


def test_model_beta_exposes_temporal_thresholds_to_modelling_backend(phase_root):
    defaults = _text(phase_root / "+phase_model_beta" / "defaultConfig.m")
    schema = _text(phase_root / "+phase_model_beta" / "schema.m")
    engine = _text(phase_root / "+phase_model_beta" / "LegacyEngine.m")
    modelling = _text(phase_root / "MatlabFunctions" / "ModellingInTime.m")

    thresholds = (
        "min_period_days",
        "min_coll_snr",
        "min_coll_corr_samples",
        "spline_min_knot_intervals",
        "spline_max_fraction",
    )
    for name in thresholds:
        assert f"'{name}_method'" in defaults
        assert f"field('{name}_method'" in schema
        assert f"'{name}'" in engine
        assert f"addParameter(p, '{name}'" in modelling
    assert "NigNoiR >= min_coll_snr" in modelling
    assert "num_obs * spline_max_fraction" in modelling


def test_model_beta_skips_unresolvable_temporal_covariance_before_fitting(phase_root):
    modelling = _text(phase_root / "MatlabFunctions" / "ModellingInTime.m")
    covariance_start = modelling.index("% -- 4) Covariance modelling")
    empirical_fit = modelling.index("f1DEmpCovEst", covariance_start)
    early_guard = modelling.index("n_obs_coll < 8 || n_epochs_coll < 8", covariance_start)

    assert early_guard < empirical_fit
    assert "The spline-only temporal result is retained." in modelling
    assert "numel(tauGrid) < 3 || numel(unique(tauGrid)) < 3" in modelling
    assert "isempty(eCovF_smooth) || isempty(tauGrid)" in modelling


def test_geosplinter_runner_closes_process_output_stream(phase_root):
    runner = _text(phase_root / "MatlabFunctions" / "runGeoSplinter.m")

    assert "process.waitFor();" in runner
    assert "reader.close();" in runner
    assert runner.index("process.waitFor();") < runner.index("reader.close();")
