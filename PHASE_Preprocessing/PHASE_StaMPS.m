function app = PHASE_StaMPS(workDir)
%PHASE_STAMPS Launch the production PHASE StaMPS application.
%
% PHASE_StaMPS() opens the selected/current ASC_* or DSC_* dataset.
% PHASE_StaMPS(WORKDIR) uses the dataset folder supplied by preprocessing.

if nargin < 1
    app = PHASE_StaMPS_beta();
else
    app = PHASE_StaMPS_beta(workDir);
end
end
