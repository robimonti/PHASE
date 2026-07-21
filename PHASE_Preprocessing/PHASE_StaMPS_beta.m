function app = PHASE_StaMPS_beta(workDir)
%PHASE_STAMPS_BETA Launch the text-based PHASE StaMPS interface.
%
%   PHASE_StaMPS_beta opens the interface for the current folder.
%   PHASE_StaMPS_beta(WORKDIR) uses an explicit ASC_*/DSC_* folder.
%
% The interface is HTML/CSS/JavaScript, while all filesystem and StaMPS
% operations remain in editable MATLAB .m files. No App Designer round-trip
% or document.xml patch is required.

launcherDir = fileparts(mfilename('fullpath'));
addpath(launcherDir);

if nargin < 1 || isempty(workDir)
    workDir = pwd;
end
workDir = char(string(workDir));

if exist(fullfile(workDir, 'input_StaMPS.mat'), 'file') ~= 2
    selected = uigetdir(workDir, ...
        'Select the ASC_*/DSC_* folder containing input_StaMPS.mat');
    if isequal(selected, 0)
        if nargout > 0
            app = [];
        end
        return
    end
    workDir = selected;
end

app = phase_stamps_beta.App(workDir, launcherDir);
if nargout == 0
    clear app
end
end
