import importlib.util
import re
import zipfile


def _text(path):
    return path.read_text(encoding="utf-8")


def _stable_xml(phase_root):
    with zipfile.ZipFile(phase_root / "legacy" / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp") as app:
        return app.read("matlab/document.xml").decode("utf-8")


def test_beta_has_a_text_launcher_and_no_new_mlapp(phase_root):
    preproc = phase_root / "PHASE_Preprocessing"
    assert (preproc / "PHASE_StaMPS_beta.m").is_file()
    assert not (preproc / "PHASE_StaMPS_beta.mlapp").exists()
    launcher = _text(preproc / "PHASE_StaMPS_beta.m")
    assert "phase_stamps_beta.App" in launcher
    assert "resolveLauncherDir" in launcher
    assert "fullfile(invocationDir,'..','PHASE_Preprocessing')" in launcher
    assert "runtimeAssetsMissing" in launcher
    assert "input_StaMPS.mat" not in launcher


def test_beta_bootstraps_missing_configuration_from_project_metadata(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta"
    detector = _text(package / "autoDetectConfig.m")
    controller = _text(package / "App.m")

    assert "projectFolder = char(java.io.File(fileparts(preprocFolder)).getCanonicalPath())" in detector
    assert "cfg.project_path = projectFolder" in detector
    assert "'project_path'" in detector
    assert "No input_StaMPS.mat was found" in controller
    assert "select the StaMPS" in controller
    assert "press Save to create it" in controller
    assert "phase_stamps_beta.findStaMPSInstallation" in detector


def test_beta_preflights_complete_stamps_runtime_before_processing(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta"
    controller = _text(package / "App.m")
    validator = _text(package / "validateConfig.m")
    inspector = _text(package / "inspectStaMPSInstallation.m")
    runtime = _text(package / "prepareRuntime.m")

    assert "phase_stamps_beta.prepareRuntime(candidate)" in controller
    assert "phase_stamps_beta.inspectStaMPSInstallation" in validator
    for executable in (
        "calamp.exe",
        "cpxsum.exe",
        "pscphase.exe",
        "pscdem.exe",
        "psclonlat.exe",
        "selpsc_patch.exe",
        "selsbc_patch.exe",
        "triangle.exe",
        "snaphu.exe",
    ):
        assert executable in inspector
    assert "addpath(genpath(stampsRoot))" in runtime
    assert "setenv('PATH'" in runtime
    assert "TRAIN runtime ready" in runtime


def test_beta_default_config_is_fully_represented_in_ui_schema(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta"
    defaults = set(re.findall(r"^cfg\.(\w+)\s*=", _text(package / "defaultConfig.m"), re.M))
    schema_fields = set(re.findall(r"item\('([^']+)'", _text(package / "schema.m")))
    synthetic = {
        "stamps_preparation": "prepare_data",
        "train_flag": "train_enabled",
        "year_0": "initial_date",
        "month_0": "initial_date",
        "day_0": "initial_date",
        "ref_centre_lonlat": "ref_centre_lon",
        "ref_centre_lonlat_w": "ref_centre_lon_w",
    }
    missing = {field for field in defaults if field not in schema_fields and field not in synthetic}
    assert not missing
    assert set(synthetic.values()) <= schema_fields
    assert "density_rand" in schema_fields


def test_beta_backend_is_reproducibly_extracted_from_stable_app(phase_root):
    tool_path = phase_root / "tools" / "extract_stamps_beta_backend.py"
    spec = importlib.util.spec_from_file_location("extract_stamps_beta_backend", tool_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    generated = _text(
        phase_root
        / "PHASE_Preprocessing"
        / "+phase_stamps_beta"
        / "runProcessing.m"
    )
    assert generated == module.extract(_stable_xml(phase_root))


def test_beta_backend_preserves_stable_setparm_and_stamps_calls(phase_root):
    stable = _stable_xml(phase_root)
    backend = _text(
        phase_root
        / "PHASE_Preprocessing"
        / "+phase_stamps_beta"
        / "runProcessing.m"
    )
    callback = stable[
        stable.index("function StartButtonPushed(app, event)") :
        stable.index("% Button pushed function: LoadButton")
    ]
    pattern = r"setparm(?:_aps)?\([^\n]+"
    backend_calls = [call.rstrip() for call in re.findall(pattern, backend)]
    stable_calls = [call.rstrip() for call in re.findall(pattern, callback)]
    assert [
        call for call in backend_calls
        if not call.startswith("setparm('density_rand'")
    ] == stable_calls
    assert "if strcmpi(select_method, 'DENSITY')" in backend
    assert "setparm('density_rand', density_rand)" in backend
    assert [call.rstrip() for call in re.findall(r"stamps\([^\n]+", backend)] == [
        call.rstrip() for call in re.findall(r"stamps\([^\n]+", callback)
    ]
    assert "stamps(stamps_first_step, stamps_last_step); % execute the complete range" in backend
    assert "stamps(stamps_first_step,7)" not in backend


def test_beta_controller_requires_save_before_start(phase_root):
    controller = _text(
        phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta" / "App.m"
    )
    assert "phase_stamps_beta.configsEqual(candidate, obj.SavedConfig)" in controller
    assert "Press Save before Start" in controller
    assert "phase_stamps_beta.runProcessing(adapter)" in controller


def test_beta_offers_model_launch_only_after_success(phase_root):
    controller = _text(
        phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta" / "App.m"
    )
    assert "if result.ok\n                obj.offerModelLaunch();" in controller
    assert "uiconfirm(obj.UIFigure" in controller
    assert "Open PHASE Model" in controller
    assert "PHASE_Model_beta();" in controller


def test_beta_exports_only_ps_time_series_not_atmospheric_delay(phase_root):
    backend = _text(
        phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta" / "runProcessing.m"
    )
    assert "Displacement time series export started" in backend
    assert "ps_plot('v-dao'" in backend
    assert "Atmospheric delay export intentionally omitted." in backend
    assert "Atmosphere time series export" not in backend
    assert "_ATMOSPHERE.xlsx" not in backend
    assert "_ATMOSPHERE.csv" not in backend
    assert "save(strcat('Atmosphere_'" not in backend


def test_beta_run_monitor_streams_matlab_and_external_output_without_cmd_windows(phase_root):
    preprocessing = phase_root / "PHASE_Preprocessing"
    package = preprocessing / "+phase_stamps_beta"
    controller = _text(package / "App.m")
    backend = _text(package / "runProcessing.m")
    runner = _text(package / "runCommandHidden.m")
    html = _text(preprocessing / "phase_stamps_beta_ui" / "index.html")
    js = _text(preprocessing / "phase_stamps_beta_ui" / "app.js")

    assert "diary(obj.LiveLogFile)" in controller
    assert "state.liveLogUrl = obj.LiveLogUrl" in controller
    assert "runtime_logs/" in controller
    assert "fetch(`${url}${separator}v=${Date.now()}`" in js
    assert "window.setInterval(pollLiveLog, 500)" in js
    assert "blocking-overlay" not in html + js
    assert "java.lang.ProcessBuilder" in runner
    assert "redirectErrorStream(true)" in runner
    assert "java.lang.String(getenv('PATH'))" in runner
    assert "pathKey = 'Path'" in runner
    assert "cmd.exe" in runner
    assert "phase_stamps_beta.runCommandHidden" in backend
    assert "mt_prep_snap_status = system" not in backend
    assert "dos('where snaphu')" not in backend
    assert "PHASE Run monitor" in backend


def test_beta_range_calendar_and_ts_picker_are_native_to_the_new_app(phase_root):
    preprocessing = phase_root / "PHASE_Preprocessing"
    package = preprocessing / "+phase_stamps_beta"
    schema = _text(package / "schema.m")
    controller = _text(package / "App.m")
    picker = _text(package / "openTsPicker.m")
    html = _text(preprocessing / "phase_stamps_beta_ui" / "index.html")
    js = _text(preprocessing / "phase_stamps_beta_ui" / "app.js")
    css = _text(preprocessing / "phase_stamps_beta_ui" / "styles.css")

    assert "item('master_date', 'Master date', 'project', 'date'" in schema
    assert 'control.type = "date"' in js
    assert "dateInputValue" in js
    assert "compactDate" in js
    assert 'style.marginLeft = `${left}%`' in js
    assert "margin-left .25s ease" in css
    assert "TSPickerOverlay" in controller
    assert "TSPickerContainer" in controller
    assert "obj.UIFigure.AutoResizeChildren = 'off'" in controller
    assert "phase_stamps_beta.openTsPicker" in controller
    assert "ts_export_picker(workDir, parentContainer" in picker
    assert "uifigure(" not in picker
    assert "No legacy MLAPP is opened" in html


def test_beta_provides_a_matlab_side_smoke_test(phase_root):
    self_test = _text(
        phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta" / "selfTest.m"
    )
    assert "phase_stamps_beta.configToUi" in self_test
    assert "phase_stamps_beta.uiToConfig" in self_test
    assert "phase_stamps_beta.saveConfig" in self_test
    assert "phase_stamps_beta.loadConfig" in self_test


def test_beta_ui_uses_local_assets_and_matlab_events(phase_root):
    ui = phase_root / "PHASE_Preprocessing" / "phase_stamps_beta_ui"
    html = _text(ui / "index.html")
    js = _text(ui / "app.js")
    css = _text(ui / "styles.css")
    assert '<link rel="stylesheet" href="styles.css">' in html
    assert '<script src="app.js"></script>' in html
    for logo in ("PHASE_logo.png", "PHASE_mod1b.png"):
        assert (ui / "assets" / logo).is_file()
        assert f'assets/{logo}' in html
    for phase_colour in ("rgb(53, 101, 207)", "rgb(203, 46, 108)", "rgb(69, 70, 70)"):
        assert phase_colour in css
    assert '"SF Pro Display"' in css
    assert "color-scheme: light" in css
    assert "http://" not in html + js + css
    assert "https://" not in html + js + css
    for event in ("Ready", "Load", "Save", "Start", "Browse", "OpenTsPicker"):
        assert f'"{event}"' in js
    assert "sendEventToMATLAB" in js
    assert 'addEventListener("DataChanged"' in js
    assert 'addEventListener("PhaseLog"' in js
    assert 'item.id === "percent_rand"' in js
    assert 'item.id === "density_rand"' in js
    assert 'item.id === "select_method"' in js
    assert '["run", "ts"].includes(PhaseUI.active)' in js
    assert "TRAIN tropospheric correction" in html
    assert '"Not applied"' in js


def test_beta_contains_every_primary_user_action(phase_root):
    html = _text(
        phase_root / "PHASE_Preprocessing" / "phase_stamps_beta_ui" / "index.html"
    )
    for element_id in (
        "load-button",
        "save-button",
        "start-button",
        "workdir",
        "open-error-log",
        "open-ts-picker",
    ):
        assert f'id="{element_id}"' in html


def test_beta_legacy_adapter_covers_all_engine_app_properties(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta"
    backend = _text(package / "runProcessing.m")
    adapter = _text(package / "LegacyAppAdapter.m")
    referenced = set(re.findall(r"app\.(\w+)", backend))
    helper_methods = {"log", "openTsPicker"}
    declared = set(re.findall(r"^\s{8}(\w+)\s*$", adapter, re.M)) | helper_methods
    assert referenced <= declared
