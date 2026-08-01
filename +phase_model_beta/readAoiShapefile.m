function [lonlat,segments,info] = readAoiShapefile(filepathAOI,filepathIN)
%READAOISHAPEFILE Read the selected AOI as geographic polygon segments.
%
% The standalone map and the numerical engine call this same function so
% the geometry displayed to the user is the geometry used to select PS.
% Projected shapefiles retain the legacy assumption that they use the UTM
% zone of the displacement dataset.

if nargin < 2, filepathIN = ''; end
filepathAOI = char(string(filepathAOI));
filepathIN = char(string(filepathIN));
if isempty(strtrim(filepathAOI)) || ~isfile(filepathAOI)
    error('PHASE_Model_beta:shapefileMissing', ...
        'The AOI shapefile does not exist: %s',filepathAOI);
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

isGeographic = all(allCoordinates(:,1)>=-180 & allCoordinates(:,1)<=180) && ...
    all(allCoordinates(:,2)>=-90 & allCoordinates(:,2)<=90);
if isGeographic
    segments = rawSegments;
    coordinateType = 'geographic';
else
    if isempty(strtrim(filepathIN)) || ~isfile(filepathIN)
        error('PHASE_Model_beta:projectedShapefileNeedsPsFile', ...
            ['The AOI shapefile uses projected coordinates. Select the displacement ', ...
             'file first so PHASE can determine its UTM zone.']);
    end
    bounds = phase_model_beta.estimatePsBoundingBox(filepathIN);
    centreLon = mean(bounds(1:2));
    centreLat = mean(bounds(3:4));
    [~,~,utmZone] = deg2utm(centreLat,centreLon);
    utmZone = utmZone(1,:);
    segments = cell(size(rawSegments));
    for partIndex = 1:numel(rawSegments)
        xy = rawSegments{partIndex};
        [lat,lon] = utm2deg(xy(:,1),xy(:,2), ...
            repmat(utmZone,size(xy,1),1));
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
