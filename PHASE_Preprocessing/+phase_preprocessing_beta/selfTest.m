function result = selfTest()
%SELFTEST Exercise the text configuration layer without running SNAP.

rootDir = tempname; mkdir(rootDir);
cleanup = onCleanup(@() rmdir(rootDir, 's')); %#ok<NASGU>
cfg = phase_preprocessing_beta.defaultConfig();
ui = phase_preprocessing_beta.configToUi(cfg);
roundTrip = phase_preprocessing_beta.uiToConfig(ui, cfg);
assert(phase_preprocessing_beta.configsEqual(cfg, roundTrip));
resumeUi = ui;
resumeUi.first_step = 3;
resumeUi.process_master = true;
resumeCfg = phase_preprocessing_beta.uiToConfig(resumeUi, cfg);
assert(~resumeCfg.process_master, ...
    'Master processing was not disabled when resuming from slave step 3.');
[resumeErrors, resumeWarnings] = ...
    phase_preprocessing_beta.validateConfig(resumeCfg, false);
assert(isempty(resumeErrors));
assert(any(contains(resumeWarnings, 'master processing is disabled')));
pathValue = phase_preprocessing_beta.saveConfig(rootDir, cfg);
assert(exist(pathValue, 'file') == 2);
[loaded, info] = phase_preprocessing_beta.loadConfig(rootDir);
assert(info.exists);
assert(phase_preprocessing_beta.configsEqual(cfg, loaded));
[errors, ~] = phase_preprocessing_beta.validateConfig(cfg, false);
assert(isempty(errors));
assert(phase_preprocessing_beta.estimateEpsg(9.8,45.0,10.2,45.4) == 32632);
assert(phase_preprocessing_beta.estimateEpsg(17.8,-34.2,18.8,-33.7) == 32734);
assert(phase_preprocessing_beta.estimateEpsg(-45,84.5,-40,86) == 3413);
selfTestGptResolver(rootDir);
selfTestReprojection(rootDir);
[pythonExecutable, pythonInfo] = phase_preprocessing_beta.resolvePython(cfg.python);
assert(startsWith(pythonInfo.version, '3.'));
result = struct('ok', true, ...
    'message', 'PHASE preprocessing beta configuration and reprojection self-test passed.', ...
    'python', pythonExecutable, 'pythonVersion', pythonInfo.version);
disp(result.message);
fprintf('Python %s: %s\n', pythonInfo.version, pythonExecutable);
end

function selfTestGptResolver(rootDir)
folder = fullfile(rootDir,'snap','bin');
mkdir(folder);
pathValue = fullfile(folder,'gpt.exe');
fileId = fopen(pathValue,'w');
assert(fileId ~= -1);
fclose(fileId);
[resolved,info] = phase_preprocessing_beta.resolveGpt(fullfile(folder,'gpt'));
assert(info.exists);
assert(endsWith(lower(resolved),'gpt.exe'));
end

function selfTestReprojection(rootDir)
pathValue = fullfile(rootDir,'reprojection_test.tif');
source = single(peaks(24));
sourceRef = georefcells([45.0 45.2],[9.8 10.2],size(source), ...
    'ColumnsStartFrom','north');
geotiffwrite(pathValue,source,sourceRef,'CoordRefSysCode',4326);
phase_preprocessing_beta.reprojectGeoTiff(pathValue,32632);
[warped,projectedRef] = readgeoraster(pathValue);
assert(isprop(projectedRef,'XWorldLimits'));
assert(~isempty(projectedRef.ProjectedCRS));
assert(isequal(projectedRef.ProjectedCRS,projcrs(32632)));
assert(any(isfinite(warped(:))));
end
