function [project, report] = importLegacy(legacyRoot, projectRoot, name, mode)
%IMPORTLEGACY Create a project and catalog final legacy results.
% MODE is 'copy' (portable result snapshot) or 'reference' (no data copy).
% Raw SAR data and partial processing state are never migrated.

if nargin < 3, name = ''; end
if nargin < 4 || isempty(mode), mode = 'copy'; end
mode = char(string(mode));
if ~any(strcmp(mode,{'copy','reference'}))
    error('PHASE:ImportModeInvalid','Mode must be copy or reference.');
end
target = phase_project.paths(projectRoot);
inventory = phase_project.scanLegacy(legacyRoot,target);
if isempty(inventory.items)
    error('PHASE:NoLegacyResults', ...
        'No Model output or StaMPS EXPORT products found in %s.',inventory.sourceRoot);
end
source = inventory.sourceRoot;
targetForComparison = target.root;
sourceForComparison = source;
if ispc
    targetForComparison = lower(targetForComparison);
    sourceForComparison = lower(sourceForComparison);
end
if strcmp(targetForComparison,sourceForComparison) || ...
        startsWith(targetForComparison,[sourceForComparison filesep])
    error('PHASE:ImportTargetInsideSource', ...
        'Choose a new project folder outside the legacy workspace.');
end
if isempty(name)
    [~,name] = fileparts(source);
end
[project,p] = phase_project.create(target.root,name);
items = inventory.items;
for k = 1:numel(items)
    items(k).status = mode;
    if strcmp(mode,'copy')
        destination = fullfile(p.root,items(k).relativeDestination);
        destinationFolder = fileparts(destination);
        if ~isfolder(destinationFolder), mkdir(destinationFolder); end
        [ok,message] = copyfile(items(k).source,destination);
        if ~ok
            items(k).status = 'copy-failed';
            warning('PHASE:LegacyCopyFailed','%s: %s',items(k).source,message);
        end
    end
end
report = struct( ...
    'schemaVersion',1, ...
    'sourceRoot',source, ...
    'mode',mode, ...
    'importedAt',char(datetime('now','TimeZone','UTC', ...
        'Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
    'modelRuns',inventory.modelRuns, ...
    'stampsDatasets',inventory.stampsDatasets, ...
    'totalBytes',inventory.totalBytes, ...
    'items',items);
phase_project.writeJson(fullfile(p.imports,'legacy-import.json'),report);
end
