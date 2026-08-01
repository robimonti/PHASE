function result = geoSplinterSelfTest()
%GEOSPLINTERSELFTEST Exercise the native analysis binary on a tiny dataset.

runtimeRoot = phase_model_beta.projectRoot();
functionsDir = fullfile(runtimeRoot,'MatlabFunctions');
addpath(functionsDir);
previousFolder = pwd;
folderRestore = onCleanup(@() cd(previousFolder)); %#ok<NASGU>
cd(runtimeRoot);

relativeRoot = fullfile('phase_model_beta_ui', ...
    ['runtime_test_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
absoluteRoot = fullfile(runtimeRoot,relativeRoot);
cleanup = onCleanup(@() removeTestFolder(absoluteRoot)); %#ok<NASGU>
inputDir = fullfile(relativeRoot,'data_input');
outputDir = fullfile(relativeRoot,'data_output');
jobDir = fullfile(relativeRoot,'job');
mkdir(fullfile(runtimeRoot,inputDir));
mkdir(fullfile(runtimeRoot,outputDir));
mkdir(fullfile(runtimeRoot,jobDir));

% Mirror the short names used by the real temporal workflow. The 2005
% Windows binary has fixed-width filename buffers for some diagnostics.
inputName = 'PS_1.txt';
outputName = 'PS_1_cub';
writematrix([0 0; 1 1; 2 0.2; 3 1.1; 4 0.1; 5 0.9], ...
    fullfile(runtimeRoot,inputDir,inputName),'Delimiter','space');
jobFile_analysis(1,2,inputName,outputName,6,2,0,5,0,10, ...
    inputDir,outputDir,jobDir);

switch detectOS()
    case 'windows'
        executable = fullfile('geoSplinter','windows','geoSplinter_analysis');
    case 'linux'
        executable = fullfile('geoSplinter','linux','geoSplinter_analysis');
    case 'macos_intel'
        executable = fullfile('geoSplinter','macos_intel','geoSplinter_analysis');
    case 'macos_apple_silicon'
        executable = fullfile('geoSplinter','macos_apple_silicon','geoSplinter_analysis');
    otherwise
        error('PHASE_Model_beta:unsupportedPlatform', ...
            'Unsupported platform for the geoSplinter self-test.');
end

jobPath = fullfile(jobDir,[outputName '.job']);
[status,output] = runGeoSplinter(executable,jobPath);
requiredSuffixes = {'.out.txt','.hdr.txt','.grd.txt', ...
    '.par.txt','.std.txt','.mat.txt'};
requiredPaths = cellfun(@(suffix) ...
    fullfile(runtimeRoot,outputDir,[outputName suffix]), ...
    requiredSuffixes,'UniformOutput',false);
missing = requiredPaths(~cellfun(@isfile,requiredPaths));
if status ~= 0 || ~isempty(missing)
    error('PHASE_Model_beta:geoSplinterSelfTestFailed', ...
        ['geoSplinter smoke test failed with status %d. Missing outputs: %s. ' ...
         'Process output: %s'], ...
        status,strjoin(missing,', '),strtrim(output));
end

result = struct('ok',true,'message', ...
    'PHASE Model geoSplinter native-runtime self-test passed.', ...
    'executable',char(java.io.File(executable).getPath()));
fprintf('%s\n',result.message);
end

function removeTestFolder(folder)
try
    if isfolder(folder), rmdir(folder,'s'); end
catch
end
end
