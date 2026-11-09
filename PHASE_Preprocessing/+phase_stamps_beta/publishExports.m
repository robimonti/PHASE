function destination = publishExports(workDir)
%PUBLISHEXPORTS Copy final StaMPS products into the PHASE 7 result area.
% The native EXPORT folder remains intact for StaMPS resume and reruns.

destination = '';
projectRoot = phase_project.findRoot(workDir);
if isempty(projectRoot), return; end
[~,p] = phase_project.open(projectRoot);
source = fullfile(workDir,'EXPORT');
if ~isfolder(source)
    error('PHASE_StaMPS_beta:exportMissing', ...
        'StaMPS completed but EXPORT is missing: %s.',source);
end
[~,datasetName] = fileparts(workDir);
destination = fullfile(p.exports,datasetName);
if ~isfolder(destination), mkdir(destination); end
files = dir(source);
files = files(~[files.isdir]);
for k = 1:numel(files)
    [ok,message] = copyfile(fullfile(source,files(k).name), ...
        fullfile(destination,files(k).name),'f');
    if ~ok
        error('PHASE_StaMPS_beta:exportPublishFailed', ...
            'Could not publish %s: %s',files(k).name,message);
    end
end
% Keep the geographic selection with the delivered results, without copying
% the full radar mask or source rasters into the user-facing export folder.
areaMetadata=fullfile(workDir,'psi_area.json');
if isfile(areaMetadata)
    [ok,message]=copyfile(areaMetadata,fullfile(destination,'psi_area.json'),'f');
    if ~ok,error('PHASE:AreaPublishFailed','Could not publish PSI area metadata: %s',message);end
end
end
