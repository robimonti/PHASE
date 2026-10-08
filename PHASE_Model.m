function app = PHASE_Model(projectRoot)
%PHASE_MODEL Launch the production PHASE geospatial modelling application.
%
% The implementation retains its phase_model_beta package name for backward
% compatibility with existing configurations and tested integrations.

if nargin < 1
    app = PHASE_Model_beta();
else
    app = PHASE_Model_beta(projectRoot);
end
end
