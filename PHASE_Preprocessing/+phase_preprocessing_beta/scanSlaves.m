function files = scanSlaves(rootDir)
%SCANSLAVES Inventory Sentinel-1 ZIP and COSMO-SkyMed HDF5 inputs.

template = struct('name','','relativePath','','date','','type','', ...
    'status','','sizeBytes',0);
files = repmat(template, 0, 1);
projectFolder = fullfile(rootDir, 'PHASE_Preprocessing');
slavesFolder = fullfile(projectFolder, 'slaves');
if ~isfolder(slavesFolder), return; end

zipFiles = dir(fullfile(slavesFolder, '**', '*.zip'));
for k = 1:numel(zipFiles)
    item = template;
    item.name = zipFiles(k).name;
    item.relativePath = relativePath(fullfile(zipFiles(k).folder, zipFiles(k).name), slavesFolder);
    rawDate = sentinelDate(zipFiles(k).name);
    item.date = displayDate(rawDate);
    item.type = 'Sentinel-1';
    item.status = ternary(~isempty(rawDate) && isProcessed(projectFolder, rawDate), ...
        'Already processed', 'Ready');
    item.sizeBytes = double(zipFiles(k).bytes);
    files(end+1) = item; %#ok<AGROW>
end

h5Files = dir(fullfile(slavesFolder, '**', '*.h5'));
for k = 1:numel(h5Files)
    item = template;
    item.name = h5Files(k).name;
    item.relativePath = relativePath(fullfile(h5Files(k).folder, h5Files(k).name), slavesFolder);
    token = regexp(h5Files(k).name, '(\d{8})', 'tokens', 'once');
    if isempty(token), rawDate = ''; else, rawDate = token{1}; end
    item.date = displayDate(rawDate);
    if startsWith(upper(h5Files(k).name), 'CSG'), item.type = 'CSG'; else, item.type = 'CSK'; end
    item.status = 'Ready';
    item.sizeBytes = double(h5Files(k).bytes);
    files(end+1) = item; %#ok<AGROW>
end

if ~isempty(files)
    [~, order] = sort(lower(string({files.name})));
    files = files(order);
end
end

function value = sentinelDate(name)
token = regexp(name, '_(\d{8})T\d{6}_', 'tokens', 'once');
if isempty(token), value = ''; else, value = token{1}; end
end

function value = displayDate(raw)
if numel(raw) == 8, value = sprintf('%s-%s-%s', raw(1:4),raw(5:6),raw(7:8));
else, value = '-'; end
end

function value = isProcessed(projectFolder, dateText)
value = ~isempty(dir(fullfile(projectFolder,'coreg',['*_' dateText '.dim']))) && ...
    ~isempty(dir(fullfile(projectFolder,'ifg',['*_' dateText '.dim'])));
end

function value = relativePath(pathValue, root)
prefix = [char(root) filesep];
if startsWith(pathValue, prefix), value = pathValue(numel(prefix)+1:end); else, value = pathValue; end
end

function value = ternary(condition, yesValue, noValue)
if condition, value = yesValue; else, value = noValue; end
end
