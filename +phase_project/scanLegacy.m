function inventory = scanLegacy(legacyRoot, destinationPaths)
%SCANLEGACY Inventory final products from a pre-project PHASE workspace.
% No file is modified. Processing intermediates are intentionally excluded.

legacyRoot = char(java.io.File(char(string(legacyRoot))).getCanonicalPath());
if ~isfolder(legacyRoot)
    error('PHASE:LegacyFolderMissing','Legacy folder does not exist: %s.',legacyRoot);
end
if nargin < 2 || isempty(destinationPaths)
    resultRoot = fullfile('results');
    modelRoot = fullfile(resultRoot,'model');
    exportsRoot = fullfile(resultRoot,'exports');
else
    modelRoot = erase(destinationPaths.model,[destinationPaths.root filesep]);
    exportsRoot = erase(destinationPaths.exports,[destinationPaths.root filesep]);
end
inventory = struct('sourceRoot',legacyRoot,'items', ...
    struct('kind',{},'source',{},'relativeDestination',{},'bytes',{}), ...
    'modelRuns',0,'stampsDatasets',0,'totalBytes',0);

runs = dir(fullfile(legacyRoot,'output_*'));
runs = runs([runs.isdir]);
for k = 1:numel(runs)
    runName = runs(k).name;
    if isempty(regexp(runName,'^output_\d+$','once')), continue; end
    runRoot = fullfile(legacyRoot,runName);
    before = numel(inventory.items);
    inventory = addFiles(inventory,runRoot,'*.xlsx', ...
        fullfile(modelRoot,runName),'model-report');
    inventory = addFiles(inventory,fullfile(runRoot,'figures'),'*', ...
        fullfile(modelRoot,runName,'figures'),'model-figure');
    inventory = addFiles(inventory,fullfile(runRoot,'files','shp'),'*', ...
        fullfile(modelRoot,runName,'gis'),'model-gis');
    inventory = addFiles(inventory,fullfile(runRoot,'files','mat'),'*', ...
        fullfile(modelRoot,runName,'mat'),'model-mat');
    if numel(inventory.items) > before
        inventory.modelRuns = inventory.modelRuns + 1;
    end
end

datasets = [dir(fullfile(legacyRoot,'ASC_*')); ...
    dir(fullfile(legacyRoot,'DSC_*')); ...
    dir(fullfile(legacyRoot,'DES_*'))];
datasets = datasets([datasets.isdir]);
for k = 1:numel(datasets)
    datasetName = datasets(k).name;
    before = numel(inventory.items);
    inventory = addFiles(inventory, ...
        fullfile(legacyRoot,datasetName,'EXPORT'),'*', ...
        fullfile(exportsRoot,datasetName),'stamps-export');
    if numel(inventory.items) > before
        inventory.stampsDatasets = inventory.stampsDatasets + 1;
    end
end
end

function inventory = addFiles(inventory,folder,pattern,destination,kind)
if ~isfolder(folder), return; end
files = dir(fullfile(folder,pattern));
files = files(~[files.isdir]);
for k = 1:numel(files)
    item = struct( ...
        'kind',kind, ...
        'source',fullfile(folder,files(k).name), ...
        'relativeDestination',fullfile(destination,files(k).name), ...
        'bytes',double(files(k).bytes));
    inventory.items(end+1) = item; %#ok<AGROW>
    inventory.totalBytes = inventory.totalBytes + item.bytes;
end
end
