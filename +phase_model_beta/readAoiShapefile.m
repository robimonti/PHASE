function [lonlat,segments,info] = readAoiShapefile(filepathAOI,filepathIN)
%READAOISHAPEFILE Read the selected AOI as geographic polygon segments.
%
% The standalone map and the numerical engine call this same function so
% the geometry displayed to the user is the geometry used to select PS.
% The source CRS is read from its .prj file; it is never inferred from
% numeric coordinate ranges or from the displacement dataset.

if nargin < 2, filepathIN = ''; end
filepathAOI = char(string(filepathAOI));
filepathIN = char(string(filepathIN));
if isempty(strtrim(filepathAOI)) || ~isfile(filepathAOI)
    error('PHASE_Model_beta:shapefileMissing', ...
        'The AOI shapefile does not exist: %s',filepathAOI);
end

sourceInfo = shapeinfo(filepathAOI);
sourceCRS = sourceInfo.CoordinateReferenceSystem;
if isempty(sourceCRS)
    error('PHASE_Model_beta:missingAoiCrs', ...
        'The AOI shapefile has no readable CRS (.prj). Define its CRS before import.');
end
if isa(sourceCRS,'projcrs')
    geographicCRS = sourceCRS.GeographicCRS;
else
    geographicCRS = sourceCRS;
end
if ~isa(geographicCRS,'geocrs') || ...
        ~contains(lower(string(geographicCRS.Name)),'wgs 84')
    error('PHASE_Model_beta:unsupportedAoiDatum', ...
        'The AOI CRS must be WGS84 geographic or WGS84-based projected CRS.');
end
features = shaperead(filepathAOI);
if isempty(features) || ~isfield(features,'X') || ~isfield(features,'Y')
    error('PHASE_Model_beta:emptyShapefile', ...
        'The AOI shapefile does not contain polygon coordinates: %s',filepathAOI);
end

rawSegments = {};
allCoordinates = zeros(0,2);
for featureIndex = 1:numel(features)
    x = double(features(featureIndex).X(:));
    y = double(features(featureIndex).Y(:));
    count = min(numel(x),numel(y));
    x = x(1:count); y = y(1:count);
    finitePair = isfinite(x) & isfinite(y);
    transitions = diff([false; finitePair; false]);
    starts = find(transitions==1);
    stops = find(transitions==-1)-1;
    for partIndex = 1:numel(starts)
        coordinates = [x(starts(partIndex):stops(partIndex)), ...
            y(starts(partIndex):stops(partIndex))];
        coordinates = removeConsecutiveDuplicates(coordinates);
        if size(unique(coordinates,'rows'),1) < 3, continue; end
        rawSegments{end+1} = coordinates; %#ok<AGROW>
        allCoordinates = [allCoordinates; coordinates]; %#ok<AGROW>
    end
end
if isempty(rawSegments)
    error('PHASE_Model_beta:emptyShapefile', ...
        'The AOI shapefile contains no valid polygon with at least three vertices.');
end

if isa(sourceCRS,'geocrs')
    segments = rawSegments;
    coordinateType = 'geographic';
else
    segments = cell(size(rawSegments));
    for partIndex = 1:numel(rawSegments)
        xy = rawSegments{partIndex};
        [lat,lon] = projinv(sourceCRS,xy(:,1),xy(:,2));
        segments{partIndex} = [lon(:),lat(:)];
    end
    coordinateType = 'projected';
end

for partIndex = 1:numel(segments)
    polygon = segments{partIndex};
    polygon = polygon(all(isfinite(polygon),2),:);
    polygon = removeConsecutiveDuplicates(polygon);
    if ~isequal(polygon(1,:),polygon(end,:))
        polygon(end+1,:) = polygon(1,:); %#ok<AGROW>
    end
    segments{partIndex} = polygon;
end

lonlat = zeros(0,2);
for partIndex = 1:numel(segments)
    if ~isempty(lonlat), lonlat(end+1,:) = [NaN NaN]; end %#ok<AGROW>
    lonlat = [lonlat; segments{partIndex}]; %#ok<AGROW>
end
finiteCoordinates = lonlat(all(isfinite(lonlat),2),:);
info = struct( ...
    'coordinateType',coordinateType, ...
    'partCount',numel(segments), ...
    'bounds',[min(finiteCoordinates(:,1)),max(finiteCoordinates(:,1)), ...
              min(finiteCoordinates(:,2)),max(finiteCoordinates(:,2))]);
end

function coordinates = removeConsecutiveDuplicates(coordinates)
if size(coordinates,1)<2, return; end
keep = [true; any(diff(coordinates,1,1)~=0,2)];
coordinates = coordinates(keep,:);
end
