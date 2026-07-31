function report = selfTest()
%SELFTEST Exercise configuration conversion and MAT persistence in MATLAB.
%
% Run before the first Windows processing test:
%   addpath('C:\path\to\PHASE\PHASE_Preprocessing')
%   phase_stamps_beta.selfTest

temporaryDir = tempname;
mkdir(temporaryDir);
cleanup = onCleanup(@() removeTemporary(temporaryDir)); %#ok<NASGU>

cfg = phase_stamps_beta.defaultConfig();
cfg.installation_folder = temporaryDir;
cfg.project_path = temporaryDir;
cfg.weed_standard_dev = 37.5;
[initial, missingInfo] = phase_stamps_beta.loadConfig(temporaryDir);
assert(~missingInfo.exists, ...
    'A new StaMPS dataset should start without input_StaMPS.mat.');
assert(phase_stamps_beta.configsEqual(initial,phase_stamps_beta.defaultConfig()), ...
    'Missing configuration did not load the canonical initial values.');
ui = phase_stamps_beta.configToUi(cfg);
roundTrip = phase_stamps_beta.uiToConfig(ui, cfg);
assert(phase_stamps_beta.configsEqual(cfg, roundTrip), ...
    'UI conversion changed one or more configuration values.');

phase_stamps_beta.saveConfig(temporaryDir, cfg);
[loaded, info] = phase_stamps_beta.loadConfig(temporaryDir);
assert(info.exists, 'Saved input_StaMPS.mat was not found.');
assert(phase_stamps_beta.configsEqual(cfg, loaded), ...
    'MAT save/load round-trip changed configuration values.');

validationCfg = cfg;
% Configuration self-tests deliberately do not require an external StaMPS
% checkout. Real Start performs the strict MATLAB/binary runtime preflight.
validationCfg.installation_folder = '';
[errors, warnings] = phase_stamps_beta.validateConfig(validationCfg, false);
assert(isempty(errors), strjoin(errors, newline));

definition = phase_stamps_beta.schema();
assert(numel(definition.items) >= 60, 'The parameter schema is incomplete.');

report = struct('ok', true, 'parameters', numel(definition.items), ...
    'warnings', {warnings}, 'matPath', info.path);
disp(report)
end

function removeTemporary(folder)
try
    if isfolder(folder), rmdir(folder, 's'); end
catch
end
end
