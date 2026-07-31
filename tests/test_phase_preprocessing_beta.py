import importlib.util
import re
import zipfile


def _text(path):
    return path.read_text(encoding="utf-8")


def _stable_xml(phase_root):
    with zipfile.ZipFile(phase_root / "PHASE_Preprocessing.mlapp") as app:
        return app.read("matlab/document.xml").decode("utf-8")


def test_preprocessing_beta_has_text_launcher_and_preserves_stable_app(phase_root):
    assert (phase_root / "PHASE_Preprocessing.mlapp").is_file()
    assert (phase_root / "PHASE_Preprocessing_beta.m").is_file()
    assert not (phase_root / "PHASE_Preprocessing_beta.mlapp").exists()
    assert "phase_preprocessing_beta.App" in _text(phase_root / "PHASE_Preprocessing_beta.m")


def test_preprocessing_engine_is_reproducibly_extracted_from_mlapp(phase_root):
    tool_path = phase_root / "tools" / "extract_preprocessing_beta_engine.py"
    spec = importlib.util.spec_from_file_location("extract_preprocessing_beta_engine", tool_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    generated = _text(
        phase_root
        / "PHASE_Preprocessing"
        / "+phase_preprocessing_beta"
        / "LegacyEngine.m"
    )
    assert generated == module.extract(_stable_xml(phase_root))


def test_preprocessing_engine_exposes_every_stable_callback_as_text(phase_root):
    stable = _stable_xml(phase_root)
    engine = _text(
        phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta" / "LegacyEngine.m"
    )
    stable_functions = re.findall(r"^\s*function\s+(?:\[[^]]+\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(", stable, re.M)
    engine_functions = re.findall(r"^\s*function\s+(?:\[[^]]+\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(", engine, re.M)
    assert len(stable_functions) == len(engine_functions)
    for callback in (
        "StartButtonPushed",
        "SearchASFButtonPushed",
        "DownloadSelectedButtonPushed",
        "ImportSENButtonPushed",
        "ImportCSKButtonPushed",
        "runSENUpdateAfterDownload",
        "initializeDownloaderMap",
    ):
        assert callback in engine_functions
    assert engine.count(
        "phase_preprocessing_beta.launchStampsBeta(stamps_app_file, stamps_app_full);"
    ) == 2
    assert "canonical_stamps_app" not in engine
    assert "dest_mlapp" not in engine
    assert "copyfile(canonical_stamps_app" not in engine


def test_preprocessing_default_config_is_fully_represented_in_schema(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    defaults = set(re.findall(r"^cfg\.(\w+)\s*=", _text(package / "defaultConfig.m"), re.M))
    schema_fields = set(re.findall(r"item\('([^']+)'", _text(package / "schema.m")))
    assert defaults == schema_fields


def test_preprocessing_beta_resolves_python3_before_running_backends(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    resolver = _text(package / "resolvePython.m")
    controller = _text(package / "App.m")
    self_test = _text(package / "selfTest.m")

    assert "%APPDATA%\\PHASE\\python.txt" in resolver
    assert "fullfile(appData, 'PHASE', 'python.txt')" in resolver
    assert "'py -3'" in resolver
    assert "major == 3 && minor >= 8" in resolver
    assert "endsWith(lower(resolved),'pythonw.exe')" in resolver
    assert "fullfile(fileparts(resolved),'python.exe')" in resolver
    assert "python3NotFound" in resolver
    assert "phase_preprocessing_beta.resolvePython(obj.Config.python)" in controller
    assert "parts = [{pythonExecutable},{scriptPath},arguments]" in controller
    assert "phase_preprocessing_beta.resolvePython(cfg.python)" in self_test


def test_preprocessing_workflow_starts_with_images_and_places_map_in_aoi(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    schema = _text(package / "schema.m")
    html = _text(
        phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui" / "index.html"
    )
    js = _text(
        phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui" / "app.js"
    )

    assert "'id', {'images','setup','aoi'" in schema
    assert "item('constellation','Satellite constellation','images'" in schema
    assert "item('update_processed_data','Update an existing processed stack','images'" in schema
    assert 'id="images-panel"' in html
    assert 'id="aoi-map-panel"' in html
    assert html.index('id="aoi-map-panel"') < html.index('id="form-grid"')
    assert 'id="data-page"' not in html
    assert "Data workspace" not in html + js
    assert 'active: "images"' in js
    assert 'group.id === "aoi"' in js
    assert 'group.id === "images"' in js
    assert 'classList.toggle("hidden", !isSEN)' in js


def test_preprocessing_controller_requires_save_and_uses_extracted_engine(phase_root):
    controller = _text(
        phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta" / "App.m"
    )
    assert "phase_preprocessing_beta.ProcessingEngine" in controller
    assert "phase_preprocessing_beta.configsEqual(candidate, obj.SavedConfig)" in controller
    assert "Press Save before Start" in controller
    assert "obj.Engine.StartButtonPushed([])" in controller
    assert "obj.Engine.StopFlag = true" in controller
    assert "obj.Engine.UIFigure.Visible = 'on'" not in controller
    assert "case 'openadvanced'" not in controller


def test_preprocessing_beta_ui_uses_local_code_phase_module_1a_and_cached_esri_tiles(phase_root):
    ui = phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui"
    html = _text(ui / "index.html")
    js = _text(ui / "app.js")
    map_js = _text(ui / "map.js")
    css = _text(ui / "styles.css")
    for asset in ("PHASE_logo.png", "PHASE_mod1a.png"):
        assert (ui / "assets" / asset).is_file()
        assert f"assets/{asset}" in html
    assert '<link rel="stylesheet" href="styles.css">' in html
    assert '<script src="map.js"></script>' in html
    assert '<script src="app.js"></script>' in html
    assert "https://" not in html + js + css
    assert "https://" not in map_js
    assert "map_tiles/" in map_js
    assert "onTilesRequested" in map_js
    assert "retryTiles(keys)" in map_js
    tile_cache = _text(phase_root / "pythonScripts" / "cache_phase_map_tiles.py")
    assert "https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile" in tile_cache
    assert "https://services.arcgisonline.com/ArcGIS/rest/services/Reference/World_Boundaries_and_Places/MapServer/tile" in tile_cache
    assert (ui / "map_tiles" / ".gitignore").is_file()
    for colour in ("rgb(53, 101, 207)", "rgb(203, 46, 108)", "rgb(69, 70, 70)"):
        assert colour in css
    assert '"SF Pro Display"' in css


def test_preprocessing_ui_contains_primary_and_advanced_actions(phase_root):
    ui = phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui"
    html = _text(ui / "index.html")
    js = _text(ui / "app.js")
    for element_id in (
        "load-button",
        "save-button",
        "start-button",
        "stop-button",
        "aoi-map",
        "map-draw",
        "map-finish",
        "map-clear",
        "map-fit",
        "map-refresh",
        "open-downloader",
        "open-update",
        "import-images",
        "refresh-slaves",
        "slaves-table-body",
        "download-map",
        "search-asf",
        "asf-results-body",
        "update-results-body",
    ):
        assert f'id="{element_id}"' in html
    for event in (
        "Ready",
        "Load",
        "Save",
        "Start",
        "Stop",
        "Browse",
        "ImportImages",
        "OpenSlavesFolder",
        "RefreshSlaves",
        "MapAoiChanged",
        "RefreshMap",
        "DownloadAoiChanged",
        "DownloadSearch",
        "DownloadLogin",
        "DownloadLogout",
        "DownloadSelected",
        "StopDownload",
        "RefreshUpdate",
        "SearchUpdate",
        "DownloadUpdate",
        "MapTilesRequested",
    ):
        assert f'"{event}"' in js
    assert "OpenAdvanced" not in js
    assert "openAdvanced" not in js


def test_preprocessing_beta_map_is_embedded_and_updates_both_aoi_contracts(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    controller = _text(package / "App.m")
    footprints = _text(package / "collectFootprints.m")
    map_base = _text(package / "mapBase.m")
    map_js = _text(
        phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui" / "map.js"
    )

    assert "phase_preprocessing_beta.mapBase()" in controller
    assert "phase_preprocessing_beta.collectFootprints" in controller
    assert "obj.Config.lon_min = values(1)" in controller
    assert "obj.Config.lon_max = values(2)" in controller
    assert 'obj.Engine.DownloaderAOIType = "polygon"' in controller
    assert "obj.Engine.DownloaderPolygonCoords = obj.MapPolygon" in controller
    assert "obj.Engine.DownloaderAOICorners = []" in controller
    assert "search_summary.json" not in footprints
    assert "downloadasf" not in footprints
    assert "fullfile(rootDir, 'PHASE_Preprocessing', 'slaves')" in footprints
    assert "manifest.safe" in footprints
    assert "Estimated Top Left Geodetic Coordinates" in footprints
    assert "load('coastlines')" in map_base
    for capability in (
        "fitToFootprints",
        "startDrawing",
        "finishDrawing",
        "vertexDrag",
        "onPolygonChanged",
        "project(lon, lat)",
        "unproject(x, y)",
        "renderBasemap",
        "localTileUrl",
        "queueTileRequest",
        "retryTiles",
    ):
        assert capability in map_js


def test_betas_have_no_runtime_mlapp_dependency(phase_root):
    preprocessing = phase_root / "PHASE_Preprocessing"
    controller = _text(preprocessing / "+phase_preprocessing_beta" / "App.m")
    engine = _text(preprocessing / "+phase_preprocessing_beta" / "LegacyEngine.m")
    processing_engine = _text(
        preprocessing / "+phase_preprocessing_beta" / "ProcessingEngine.m"
    )
    stamps_launcher = _text(preprocessing / "PHASE_StaMPS_beta.m")
    stamps_controller = _text(preprocessing / "+phase_stamps_beta" / "App.m")
    stamps_backend = _text(preprocessing / "+phase_stamps_beta" / "runProcessing.m")

    assert "phase_preprocessing_beta.ProcessingEngine" in controller
    assert "phase_preprocessing_beta.LegacyEngine" in processing_engine
    assert "PHASE_Preprocessing.mlapp" not in controller
    assert "canonical_stamps_app" not in engine
    assert "dest_mlapp" not in engine
    assert "phase_preprocessing_beta.launchStampsBeta" in engine
    assert "phase_stamps_beta.App" in stamps_launcher
    assert "phase_stamps_beta.runProcessing" in stamps_controller
    assert "run(" not in stamps_launcher
    assert "PHASE_StaMPS.mlapp" not in stamps_launcher
    assert "run(stamps_app_file)" not in stamps_backend


def test_preprocessing_dependencies_step_range_cleanup_and_output_contract(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    schema = _text(package / "schema.m")
    defaults = _text(package / "defaultConfig.m")
    controller = _text(package / "App.m")
    ui = phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui"
    js = _text(ui / "app.js")
    css = _text(ui / "styles.css")

    assert 'setFieldDisabled("master_date", Boolean(PhaseUI.config.auto_master))' in js
    assert 'item.id === "first_step"' in js
    assert '["first_step", "num_gcp"].includes(item.id)' in js
    assert 'const resumeSlaves = Number(PhaseUI.config.first_step || 1) > 1' in js
    assert "PhaseUI.config.process_master = false" in js
    assert 'setFieldDisabled("process_master", resumeSlaves)' in js
    assert 'id="resume-warning"' in _text(ui / "index.html")
    assert '.field[data-field="remove_slaves_after_processing"] { grid-column: 1; }' in css
    assert "style.marginLeft" in js
    for description in (
        "Prepare slave acquisitions",
        "Coregister and form interferograms",
        "Export the stack to StaMPS",
        "Terrain-correct optional products",
    ):
        assert description in schema
    for field in (
        "remove_slaves_after_processing",
        "remove_split_after_processing",
        "remove_coreg_after_processing",
        "remove_ifg_after_processing",
        "auto_epsg",
        "generate_coherence",
        "generate_lia",
    ):
        assert f"cfg.{field}" in defaults
        assert f"item('{field}'" in schema
    assert "original ZIP/HDF5 files" in schema
    assert "finalizeBetaProducts" in controller
    assert "removeSourceImages" in controller
    assert "if obj.StopRequested" in controller
    assert "cleanup was skipped" in controller
    assert "removeProcessingFolder('coreg'" in controller
    assert "removeProcessingFolder('ifg'" in controller
    apply_config = _text(package / "applyConfig.m")
    assert "engine.slaves_removal_SEN = 1" in apply_config
    assert "engine.slaves_removal_CSK = 1" in apply_config
    assert ".field.inactive" in css


def test_preprocessing_live_monitor_and_silent_process_contract(phase_root):
    preprocessing = phase_root / "PHASE_Preprocessing"
    package = preprocessing / "+phase_preprocessing_beta"
    controller = _text(package / "App.m")
    engine = _text(package / "LegacyEngine.m")
    runner = _text(package / "runCommandLive.m")
    resolver = _text(package / "resolveGpt.m")
    schema = _text(package / "schema.m")
    html = _text(preprocessing / "phase_preprocessing_beta_ui" / "index.html")
    js = _text(preprocessing / "phase_preprocessing_beta_ui" / "app.js")
    css = _text(preprocessing / "phase_preprocessing_beta_ui" / "styles.css")
    python_helper = _text(preprocessing / "snap2stamps" / "bin" / "phase_subprocess.py")

    assert "ExternalProgressCallback" in engine
    assert "phase_preprocessing_beta.runCommandLive" in engine
    assert "snap2stamps_slaves.bat &" not in engine
    assert "open -a Terminal ' path_2_slaves" not in engine
    assert "xterm space path_2_slaves" not in engine
    assert "java.lang.ProcessBuilder" in runner
    assert "redirectErrorStream(true)" in runner
    assert "PYTHONUNBUFFERED" in runner
    assert "engine.updateOutput(line)" in runner
    assert "PhaseRunProgress" in controller + js
    assert "fprintf('[PHASE preprocessing beta" not in controller
    assert "processing-progress-bar" in html
    assert "processing-progress-card" in css
    assert "blocking-overlay" not in html + js
    assert "advanced-toggle" not in html + js
    assert "'resources'" not in schema
    assert "item('cpu','Processing cores','setup'" in schema
    assert "[configured '.exe']" in resolver
    assert "CREATE_NO_WINDOW" in python_helper
    assert "stream.write(chunk)" in python_helper
    stopper = _text(package / "forceStopProcess.m")
    processing_engine = _text(package / "ProcessingEngine.m")
    assert "ActiveProcess" in processing_engine
    assert "forceStopActiveProcess" in processing_engine
    assert "taskkill.exe" in stopper
    assert "'/T','/F'" in stopper
    assert "Force stop now" in html
    assert "forceStopActiveProcess" in controller
    assert "PHASE:ProcessingStopped" in runner


def test_preprocessing_beta_creates_dataset_and_bootstraps_unconfigured_stamps(phase_root):
    engine = _text(
        phase_root
        / "PHASE_Preprocessing"
        / "+phase_preprocessing_beta"
        / "LegacyEngine.m"
    )
    assert engine.count("if ~isfolder(stamps_folder_full)") == 2
    assert engine.count("mkdir(stamps_folder_full);") == 2
    assert engine.count("if isfile(dst_input_mat)") >= 4
    assert "StaMPS launch deferred" not in engine
    assert engine.count("choice = 'Open now';") == 2
    assert engine.count("Opening PHASE_StaMPS_beta to create the initial configuration") == 2
    assert engine.count("phase_preprocessing_beta.launchStampsBeta") == 2
    launcher = _text(
        phase_root
        / "PHASE_Preprocessing"
        / "+phase_preprocessing_beta"
        / "launchStampsBeta.m"
    )
    assert "cd(launcherDir)" in launcher
    assert "PHASE_StaMPS_beta(workDir)" in launcher


def test_preprocessing_resume_skips_master_for_sentinel_and_csk(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    ui_to_config = _text(package / "uiToConfig.m")
    load_config = _text(package / "loadConfig.m")
    validate_config = _text(package / "validateConfig.m")
    apply_config = _text(package / "applyConfig.m")

    assert "if cfg.first_step > 1" in ui_to_config
    assert "cfg.process_master = false;" in ui_to_config
    assert "if cfg.first_step > 1" in load_config
    assert "cfg.process_master = false;" in load_config
    assert "Resuming from slave step %d" in validate_config
    assert "engine.master_processing_SEN = double(~cfg.process_master)" in apply_config
    assert "engine.master_processing_CSK = double(~cfg.process_master)" in apply_config


def test_preprocessing_uses_true_raster_warp_and_aoi_crs_estimator(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    engine = _text(package / "ProcessingEngine.m")
    warp = _text(package / "reprojectGeoTiff.m")
    estimator = _text(package / "estimateEpsg.m")
    controller = _text(package / "App.m")
    js = _text(
        phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui" / "app.js"
    )

    assert "classdef ProcessingEngine < phase_preprocessing_beta.LegacyEngine" in engine
    assert "phase_preprocessing_beta.reprojectGeoTiff" in engine
    for operation in (
        "readgeoraster",
        "projfwd",
        "projinv",
        "worldGrid",
        "geointerp",
        "mapinterp",
        "atomicGeoTiffWrite",
    ):
        assert operation in warp
    assert "geotiffinfo" not in warp
    assert "geotiffread" not in warp
    assert "blockRows" in warp
    assert "phase_preprocessing_beta.estimateEpsg" in controller
    assert "phase_preprocessing_beta.estimateEpsg" in _text(package / "uiToConfig.m")
    assert "32600 + zone" in estimator
    assert "32700 + zone" in estimator
    assert "code = 3413" in estimator
    assert "code = 3031" in estimator
    assert "function estimateEpsg" in js


def test_preprocessing_downloader_and_import_are_integrated_without_legacy_window(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    controller = _text(package / "App.m")
    html = _text(
        phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui" / "index.html"
    )
    js = _text(
        phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui" / "app.js"
    )

    for helper in (
        "defaultDownloader.m",
        "readAsfSearch.m",
        "scanSlaves.m",
        "sentinelUpdateContext.m",
    ):
        assert (package / helper).is_file()

    assert "obj.importImages()" in controller
    assert "uigetfile({'*.zip'" in controller
    assert "uigetfile({'*.h5'" in controller
    assert "phase_preprocessing_beta.scanSlaves" in controller
    assert "controller.py" in controller
    assert "{'search'}" in controller
    assert "{'login'}" in controller
    assert "phase_download_manager.py" in controller
    assert "search_update_sentinel1_images.py" in controller
    assert "obj.Engine.UIFigure.Visible = 'on'" not in controller
    assert 'id="downloader-panel"' in html
    assert 'id="update-panel"' in html
    assert 'id="earthdata-username"' in html
    assert 'id="earthdata-password"' in html
    assert 'id="filter-processing-level"' in html
    assert 'id="filter-beam-mode"' in html
    assert 'id="filter-polarization"' in html
    assert 'id="filter-flight-direction"' in html
    assert 'id="filter-subtype"' in html
    assert 'new PhaseMap(byId("download-map")' in js
    assert 'send("DownloadSearch"' in js
    assert 'send("DownloadSelected"' in js
    assert 'send("ImportImages"' in js


def test_preprocessing_downloads_are_async_resumable_and_visible(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    controller = _text(package / "App.m")
    ui = phase_root / "PHASE_Preprocessing" / "phase_preprocessing_beta_ui"
    html = _text(ui / "index.html")
    js = _text(ui / "app.js")
    css = _text(ui / "styles.css")
    manager = _text(phase_root / "pythonScripts" / "phase_download_manager.py")

    assert "java.lang.ProcessBuilder" in controller
    assert "timer('ExecutionMode','fixedSpacing'" in controller
    assert "pollDownloadTransfer" in controller
    assert "destroyForcibly" in controller
    assert "'transfer',obj.Transfer" in controller
    assert "phase_download_manager.py" in controller
    assert "downloader_update_images.py" not in controller
    for element_id in (
        "asf-sort-key",
        "asf-sort-direction",
        "download-transfer-card",
        "download-transfer-percent",
        "download-transfer-progress",
        "download-transfer-count",
        "stop-download-transfer",
        "update-sort-key",
        "update-sort-direction",
        "update-transfer-card",
        "update-transfer-percent",
        "update-transfer-progress",
        "update-transfer-count",
        "stop-update-transfer",
    ):
        assert f'id="{element_id}"' in html
    assert "sortAsfResults" in js
    assert "sortUpdateResults" in js
    assert "renderTransferCard" in js
    assert "Force stop download" in html
    assert ".download-map-column { display: flex; flex-direction: column; }" in css
    assert "flex: 1 1 560px" in css
    assert ".transfer-percentage" in css
    assert "os.replace(partial, target)" in manager
    assert '"Range": f"bytes={offset}-"' in manager
    assert "stop_if_requested" in manager
    assert "MAX_ATTEMPTS = 3" in manager


def test_preprocessing_beta_has_matlab_side_smoke_test(phase_root):
    self_test = _text(
        phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta" / "selfTest.m"
    )
    assert "phase_preprocessing_beta.configToUi" in self_test
    assert "phase_preprocessing_beta.uiToConfig" in self_test
    assert "phase_preprocessing_beta.saveConfig" in self_test
    assert "phase_preprocessing_beta.loadConfig" in self_test
    assert "phase_preprocessing_beta.reprojectGeoTiff" in self_test
