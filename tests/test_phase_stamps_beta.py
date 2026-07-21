import importlib.util
import re
import zipfile


def _text(path):
    return path.read_text(encoding="utf-8")


def _stable_xml(phase_root):
    with zipfile.ZipFile(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp") as app:
        return app.read("matlab/document.xml").decode("utf-8")


def test_beta_has_a_text_launcher_and_no_new_mlapp(phase_root):
    preproc = phase_root / "PHASE_Preprocessing"
    assert (preproc / "PHASE_StaMPS_beta.m").is_file()
    assert not (preproc / "PHASE_StaMPS_beta.mlapp").exists()
    assert "phase_stamps_beta.App" in _text(preproc / "PHASE_StaMPS_beta.m")


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
    assert [call.rstrip() for call in re.findall(pattern, backend)] == [
        call.rstrip() for call in re.findall(pattern, callback)
    ]
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
    assert "http://" not in html + js + css
    assert "https://" not in html + js + css
    for event in ("Ready", "Load", "Save", "Start", "Browse", "OpenTsPicker"):
        assert f'"{event}"' in js
    assert "sendEventToMATLAB" in js
    assert 'addEventListener("DataChanged"' in js
    assert 'addEventListener("PhaseLog"' in js


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
