function result = selfTest()
%SELFTEST Verify that the standalone Model beta runtime is complete.

rootDir = phase_model_beta.projectRoot();
required = {
    fullfile(rootDir, 'PHASE_Model_beta.m')
    fullfile(rootDir, '+phase_model_beta', 'App.m')
    fullfile(rootDir, '+phase_model_beta', 'LegacyEngine.m')
    fullfile(rootDir, '+phase_model_beta', 'applyConfigToEngine.m')
    fullfile(rootDir, '+phase_model_beta', 'estimatePsBoundingBox.m')
    fullfile(rootDir, '+phase_model_beta', 'readAoiShapefile.m')
    fullfile(rootDir, '+phase_model_beta', 'exportFigure.m')
    fullfile(rootDir, '+phase_model_beta', 'throwIfStopped.m')
    fullfile(rootDir, '+phase_model_beta', 'mapBase.m')
    fullfile(rootDir, '+phase_model_beta', 'resolvePythonPath.m')
    fullfile(rootDir, '+phase_model_beta', 'geoSplinterSelfTest.m')
    fullfile(rootDir, 'MatlabFunctions', 'runGeoSplinter.m')
    fullfile(rootDir, 'MatlabFunctions', 'interp1Unique.m')
    fullfile(rootDir, 'geoSplinter', 'windows', 'geoSplinter_analysis.exe')
    fullfile(rootDir, 'geoSplinter', 'windows', 'geoSplinter_synthesis.exe')
    fullfile(rootDir, 'MatlabFunctions')
    fullfile(rootDir, 'PHASE_logo.png')
    fullfile(rootDir, 'phase_model_beta_ui', 'index.html')
    fullfile(rootDir, 'phase_model_beta_ui', 'app.js')
    fullfile(rootDir, 'phase_model_beta_ui', 'styles.css')
    fullfile(rootDir, 'PHASE_Preprocessing', 'phase_preprocessing_beta_ui', 'map.js')
};

missing = required(~cellfun(@(path) exist(path, 'file') == 2 || ...
    exist(path, 'dir') == 7, required));
if ~isempty(missing)
    error('PHASE_Model_beta:selfTestMissingRuntime', ...
        'Standalone Model runtime is incomplete. Missing: %s', ...
        strjoin(missing, ', '));
end

launcher = fileread(fullfile(rootDir, 'PHASE_Model_beta.m'));
engine = fileread(fullfile(rootDir, '+phase_model_beta', 'LegacyEngine.m'));
if ~contains(launcher, 'phase_model_beta.App') || ...
        ~contains(engine, 'function StartButtonPushed(app')
    error('PHASE_Model_beta:selfTestInvalidRuntime', ...
        'Standalone Model launcher or processing backend is invalid.');
end

temporaryDir = tempname;
mkdir(temporaryDir);
temporaryCleanup = onCleanup(@() removeTemporary(temporaryDir)); %#ok<NASGU>
config = phase_model_beta.defaultConfig();
config.min_period_days_method = 'manual';
config.min_period_days = 400;
config.flag_AOIbb = true;
config.aoi_polygon_lonlat = [9 45; 9.2 45.05; 9.1 45.2; 9 45];
config.lonMinAOI = 9; config.lonMaxAOI = 9.2;
config.latMinAOI = 45; config.latMaxAOI = 45.2;
phase_model_beta.saveConfig(temporaryDir,config);
[loaded,info] = phase_model_beta.loadConfig(temporaryDir);
assert(info.exists && phase_model_beta.configsEqual(config,loaded), ...
    'PHASE Model configuration save/load round-trip failed.');
[errors,~] = phase_model_beta.validateConfig(loaded,false);
assert(isempty(errors),'Standalone polygon AOI validation failed.');

samplePath = fullfile(temporaryDir,'sample_ps.xlsx');
sample = [ ...
    NaN NaN NaN NaN 0
    NaN NaN NaN NaN 1
    1 9.0 45.0 0 0
    2 11.0 47.0 0 0];
writematrix(sample,samplePath);
bounds = phase_model_beta.estimatePsBoundingBox(samplePath);
assert(bounds(1) < 9 && bounds(2) > 11 && ...
    bounds(3) < 45 && bounds(4) > 47, ...
    'Full PS extent estimation failed.');

result = struct( ...
    'ok', true, ...
    'message', 'PHASE Model beta standalone self-test passed.', ...
    'engine', fullfile(rootDir, '+phase_model_beta', 'LegacyEngine.m'));
fprintf('%s\n', result.message);
end

function removeTemporary(folder)
try
    if isfolder(folder), rmdir(folder,'s'); end
catch
end
end
