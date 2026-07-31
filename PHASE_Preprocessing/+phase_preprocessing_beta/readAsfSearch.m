function state = readAsfSearch(rootDir)
%READASFSEARCH Convert downloader JSON results to a uihtml-safe state.

template = struct('id','','sceneName','','date','','time','', ...
    'pathNumber',0,'frameNumber',0,'direction','','sizeGB',0, ...
    'sizeBytes',0,'selected',false,'coordinates',zeros(0,2));
state = struct('results',{repmat(template,0,1)},'count',0, ...
    'totalSizeGB',0,'recommended','');
summaryPath = fullfile(rootDir, 'downloadasf', 'search_summary.json');
if ~isfile(summaryPath), return; end
summary = jsondecode(fileread(summaryPath));
if isfield(summary,'status') && strcmpi(char(string(summary.status)),'invalid')
    if isfield(summary,'message'), error('PHASE:ASFInvalidResults','%s',char(string(summary.message))); end
end
if ~isfield(summary,'products') || isempty(summary.products), return; end

products = summary.products;
results = repmat(template, numel(products), 1);
for k = 1:numel(products)
    product = products(k);
    scene = char(string(product.sceneName));
    [dateText,timeText] = sceneDateTime(scene, product);
    results(k).id = ['asf-' num2str(k)];
    results(k).sceneName = scene;
    results(k).date = dateText;
    results(k).time = timeText;
    results(k).pathNumber = double(product.pathNumber);
    results(k).frameNumber = double(product.frameNumber);
    results(k).direction = char(string(product.flightDirection));
    results(k).sizeGB = double(product.size);
    results(k).sizeBytes = double(product.size_bytes);
    results(k).selected = false;
    results(k).coordinates = footprintCoordinates(product.footprint.coordinates);
end
state.results = results;
state.count = numel(results);
if isfield(summary,'total_size_gb'), state.totalSizeGB = double(summary.total_size_gb); end

paths = unique([products.pathNumber]); frames = unique([products.frameNumber]);
directions = unique(string({products.flightDirection}));
if numel(paths)==1 && numel(frames)==1 && numel(directions)==1
    state.recommended = sprintf('Compatible stack: path %d · frame %d · %s', ...
        paths(1),frames(1),char(directions(1)));
elseif isfield(summary,'best_path') && ~isempty(summary.best_path) && ...
        isfield(summary,'best_frame') && ~isempty(summary.best_frame) && ...
        isfield(summary,'best_direction') && ~isempty(summary.best_direction)
    state.recommended = sprintf('Recommended: path %d · frame %d · %s', ...
        summary.best_path,summary.best_frame,char(string(summary.best_direction)));
else
    state.recommended = 'Select products with one common path, frame and direction.';
end
end

function [dateText,timeText] = sceneDateTime(scene, product)
token = regexp(scene, '(\d{8})T(\d{6})', 'tokens', 'once');
if ~isempty(token)
    raw = token{1}; time = token{2};
    dateText = sprintf('%s-%s-%s',raw(1:4),raw(5:6),raw(7:8));
    timeText = sprintf('%s:%s:%s',time(1:2),time(3:4),time(5:6));
elseif isfield(product,'startTime')
    dateText = char(string(product.startTime)); timeText = '';
else
    dateText = ''; timeText = '';
end
end

function coords = footprintCoordinates(raw)
while iscell(raw) && ~isempty(raw), raw = raw{1}; end
coords = squeeze(double(raw));
if size(coords,2) ~= 2 && size(coords,1)==2, coords = coords.'; end
if size(coords,2) ~= 2, coords = zeros(0,2); end
end
