function applyConfigToEngine(engine,config)
%APPLYCONFIGTOENGINE Apply typed configuration without using hidden widgets.

if isempty(engine) || ~isvalid(engine)
    error('PHASE_Model_beta:engineUnavailable', ...
        'The PHASE Model processing engine is not available.');
end
if ~isstruct(config) || ~isscalar(config)
    error('PHASE_Model_beta:invalidEngineConfig', ...
        'The PHASE Model configuration must be a scalar structure.');
end

names = fieldnames(config);
for k = 1:numel(names)
    name = names{k};
    if isprop(engine,name)
        engine.(name) = config.(name);
    end
end
engine.pythonPath = phase_model_beta.resolvePythonPath(engine.pythonPath);
end
