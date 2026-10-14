function applyConfigToEngine(engine,config)
%APPLYCONFIGTOENGINE Apply typed configuration to the native model runner.

if isempty(engine) || ~isvalid(engine)
    error('PHASE_Model_beta:engineUnavailable', ...
        'The PHASE Model processing engine is not available.');
end
if ~isstruct(config) || ~isscalar(config)
    error('PHASE_Model_beta:invalidEngineConfig', ...
        'The PHASE Model configuration must be a scalar structure.');
end

engine.applyConfig(config);
end
