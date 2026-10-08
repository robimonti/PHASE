function results = importedResults(projectRoot)
%IMPORTEDRESULTS Resolve cataloged legacy results and detect missing files.

[~,p] = phase_project.open(projectRoot);
catalogPath = fullfile(p.imports,'legacy-import.json');
results = struct('kind',{},'path',{},'exists',{},'copied',{},'source',{});
if ~isfile(catalogPath), return; end
catalog = jsondecode(fileread(catalogPath));
if ~isfield(catalog,'schemaVersion') || catalog.schemaVersion ~= 1 || ...
        ~isfield(catalog,'items')
    error('PHASE:ImportCatalogInvalid','Invalid import catalog: %s.',catalogPath);
end
for k = 1:numel(catalog.items)
    item = catalog.items(k);
    copied = strcmp(char(string(item.status)),'copy');
    if copied
        pathValue = fullfile(p.root,char(string(item.relativeDestination)));
    else
        pathValue = char(string(item.source));
    end
    results(end+1) = struct( ...
        'kind',char(string(item.kind)), ...
        'path',pathValue, ...
        'exists',isfile(pathValue), ...
        'copied',copied, ...
        'source',char(string(item.source))); %#ok<AGROW>
end
end
