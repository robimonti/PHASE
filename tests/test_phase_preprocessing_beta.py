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
    assert engine.count("PHASE_StaMPS_beta(stamps_app_full);") == 2
    assert "canonical_stamps_app" not in engine
    assert "dest_mlapp" not in engine
    assert "copyfile(canonical_stamps_app" not in engine


def test_preprocessing_default_config_is_fully_represented_in_schema(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta"
    defaults = set(re.findall(r"^cfg\.(\w+)\s*=", _text(package / "defaultConfig.m"), re.M))
    schema_fields = set(re.findall(r"item\('([^']+)'", _text(package / "schema.m")))
    assert defaults == schema_fields


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
    assert "phase_preprocessing_beta.LegacyEngine" in controller
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
    stamps_launcher = _text(preprocessing / "PHASE_StaMPS_beta.m")
    stamps_controller = _text(preprocessing / "+phase_stamps_beta" / "App.m")
    stamps_backend = _text(preprocessing / "+phase_stamps_beta" / "runProcessing.m")

    assert "phase_preprocessing_beta.LegacyEngine" in controller
    assert "PHASE_Preprocessing.mlapp" not in controller
    assert "canonical_stamps_app" not in engine
    assert "dest_mlapp" not in engine
    assert "PHASE_StaMPS_beta(stamps_app_full);" in engine
    assert "phase_stamps_beta.App" in stamps_launcher
    assert "phase_stamps_beta.runProcessing" in stamps_controller
    assert "run(" not in stamps_launcher
    assert "PHASE_StaMPS.mlapp" not in stamps_launcher
    assert "run(stamps_app_file)" not in stamps_backend


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
    assert "{'download'}" in controller
    assert "search_update_sentinel1_images.py" in controller
    assert "downloader_update_images.py" in controller
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


def test_preprocessing_beta_has_matlab_side_smoke_test(phase_root):
    self_test = _text(
        phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta" / "selfTest.m"
    )
    assert "phase_preprocessing_beta.configToUi" in self_test
    assert "phase_preprocessing_beta.uiToConfig" in self_test
    assert "phase_preprocessing_beta.saveConfig" in self_test
    assert "phase_preprocessing_beta.loadConfig" in self_test
