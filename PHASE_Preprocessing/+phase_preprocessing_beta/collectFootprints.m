function footprints = collectFootprints(rootDir)
%COLLECTFOOTPRINTS Read one valid local SEN and CSK footprint from slaves.

template = struct('id', '', 'name', '', 'source', '', ...
    'selected', false, 'coordinates', zeros(0, 2));
footprints = repmat(template, 0, 1);

slavesFolder = fullfile(rootDir, 'PHASE_Preprocessing', 'slaves');
senFiles = dir(fullfile(slavesFolder, '**', '*.zip'));
for k = 1:numel(senFiles)
    try
        coords = sentinelFootprint(fullfile(senFiles(k).folder, senFiles(k).name));
        if size(coords, 1) < 3, continue; end
        item = template; item.id = ['sen-' num2str(k)];
        item.name = senFiles(k).name; item.source = 'Local Sentinel-1';
        item.selected = true; item.coordinates = coords;
        footprints(end+1) = item; %#ok<AGROW>
        break; % Match the stable map: one valid local product defines the stack footprint.
    catch
    end
end

cskFiles = dir(fullfile(slavesFolder, '**', '*.h5'));
for k = 1:numel(cskFiles)
    try
        pathValue = fullfile(cskFiles(k).folder, cskFiles(k).name);
        tl = h5readatt(pathValue, '/', 'Estimated Top Left Geodetic Coordinates');
        tr = h5readatt(pathValue, '/', 'Estimated Top Right Geodetic Coordinates');
        br = h5readatt(pathValue, '/', 'Estimated Bottom Right Geodetic Coordinates');
        bl = h5readatt(pathValue, '/', 'Estimated Bottom Left Geodetic Coordinates');
        coords = [tl(2) tl(1); tr(2) tr(1); br(2) br(1); bl(2) bl(1); tl(2) tl(1)];
        item = template; item.id = ['csk-' num2str(k)];
        item.name = cskFiles(k).name; item.source = 'Local COSMO-SkyMed';
        item.selected = true; item.coordinates = double(coords);
        footprints(end+1) = item; %#ok<AGROW>
        break; % Match the stable map: one valid local product defines the stack footprint.
    catch
    end
end
end

function coords = sentinelFootprint(pathValue)
coords = zeros(0, 2);
zipObj = java.util.zip.ZipFile(pathValue);
cleanup = onCleanup(@() zipObj.close()); %#ok<NASGU>
entries = zipObj.entries(); manifest = '';
while entries.hasMoreElements()
    entry = entries.nextElement();
    if endsWith(char(entry.getName()), 'manifest.safe')
        stream = zipObj.getInputStream(entry);
        scanner = java.util.Scanner(stream).useDelimiter('\A');
        scannerCleanup = onCleanup(@() scanner.close()); %#ok<NASGU>
        if scanner.hasNext(), manifest = char(scanner.next()); end
        break
    end
end
token = regexp(manifest, '<gml:coordinates>([^<]+)</gml:coordinates>', 'tokens', 'once');
if isempty(token), return; end
pairs = strsplit(strtrim(token{1}));
coords = zeros(numel(pairs), 2);
for k = 1:numel(pairs)
    latLon = strsplit(pairs{k}, ',');
    coords(k,:) = [str2double(latLon{2}), str2double(latLon{1})];
end
if ~isequal(coords(1,:), coords(end,:)), coords(end+1,:) = coords(1,:); end
end
