function value = displayVersion(installRoot)
%DISPLAYVERSION One version label for the hub and every workflow module.
if nargin < 1 || isempty(installRoot)
    installRoot = phase_project.installationRoot();
end
installRoot = char(string(installRoot));
value = 'v7.0.3 dev';

% Update packages carry their release tag with the engine itself.
releaseFile = fullfile(installRoot,'phase-release.json');
if isfile(releaseFile)
    try
        release = jsondecode(fileread(releaseFile));
        if isfield(release,'tag')
            value = releaseLabel(release.tag,value);
        end
    catch
    end
end

% Managed installations record the version applied by the updater.
[~,folderName] = fileparts(installRoot);
if strcmp(folderName,'engine')
    metadata = fullfile(fileparts(installRoot),'install.json');
    if isfile(metadata)
        try
            installed = jsondecode(fileread(metadata));
            if isfield(installed,'version')
                value = releaseLabel(installed.version,value);
            end
        catch
        end
    end
end
end

function value = releaseLabel(candidate,fallback)
value = fallback;
candidate = char(string(candidate));
if ~isempty(regexp(candidate,'^v?\d+\.\d+\.\d+$','once'))
    if ~startsWith(candidate,'v'), candidate = ['v' candidate]; end
    value = candidate;
end
end
