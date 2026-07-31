function [config, info] = loadConfig(rootDir)
%LOADCONFIG Load input_model.mat while filling fields added by newer betas.

if nargin < 1 || isempty(rootDir)
    rootDir = phase_model_beta.projectRoot();
end
pathValue = fullfile(rootDir,'input_model.mat');
config = phase_model_beta.defaultConfig();
info = struct('exists',false,'path',pathValue);
if ~isfile(pathValue)
    return
end

loaded = load(pathValue,'config');
if ~isfield(loaded,'config') || ~isstruct(loaded.config)
    error('PHASE_Model_beta:invalidConfig', ...
        '%s does not contain the expected config structure.',pathValue);
end
names = fieldnames(config);
for k = 1:numel(names)
    if isfield(loaded.config,names{k})
        config.(names{k}) = loaded.config.(names{k});
    end
end

legacyCovariancePairs = {
    'dtCov_method_STC1D','dtCov_STC1D'
    'dsCov_method_STC1D','dsCov_STC1D'
    'dtCov_method_STC2D','dtCov_STC2D'
    'dsCov_method_STC2D','dsCov_STC2D'
};
for k = 1:size(legacyCovariancePairs,1)
    methodName = legacyCovariancePairs{k,1};
    valueName = legacyCovariancePairs{k,2};
    legacyValue = [];
    if isfield(loaded.config,valueName)
        legacyValue = loaded.config.(valueName);
    end
    if ~isfield(loaded.config,methodName) && ...
            isscalar(legacyValue) && isfinite(legacyValue)
        config.(methodName) = 'manual';
    end
end
info.exists = true;
end
