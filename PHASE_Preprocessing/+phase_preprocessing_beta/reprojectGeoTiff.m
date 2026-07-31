function reprojectGeoTiff(filePath, epsgCode, bandIndex)
%REPROJECTGEOTIFF Warp a geographic/projected GeoTIFF into a projected CRS.
%
% Target pixel centres are inverse-projected and sampled from the source
% raster. Processing is blocked by rows to keep large products memory-safe.

if nargin < 3, bandIndex = []; end
[source, sourceRef] = readgeoraster(filePath);
if isempty(sourceRef)
    error('PHASE:MissingRasterReference', ...
        'GeoTIFF has no usable spatial reference: %s',filePath);
end
if ~isempty(bandIndex)
    if ndims(source) < 3 || size(source,3) < bandIndex
        error('PHASE:MissingRasterBand', ...
            'Band %d is missing from %s.',bandIndex,filePath);
    end
    source = source(:,:,bandIndex);
end

targetCrs = projcrs(double(epsgCode));
isGeographic = isprop(sourceRef,'LatitudeLimits') && ...
    isprop(sourceRef,'LongitudeLimits');
isProjected = isprop(sourceRef,'XWorldLimits') && ...
    isprop(sourceRef,'YWorldLimits');
if ~(isGeographic || isProjected)
    error('PHASE:UnsupportedRasterReference', ...
        'Unsupported GeoTIFF reference %s for %s.',class(sourceRef),filePath);
end

sourceCrs = [];
if isProjected
    sourceCrs = sourceRef.ProjectedCRS;
    if isempty(sourceCrs)
        error('PHASE:MissingProjectedCRS', ...
            'Projected GeoTIFF has no CRS definition: %s',filePath);
    end
    if isequal(sourceCrs,targetCrs)
        if ~isempty(bandIndex)
            atomicGeoTiffWrite(filePath,source,sourceRef,epsgCode);
        end
        return
    end
end

[boundaryLat,boundaryLon] = rasterBoundary(sourceRef,sourceCrs,isGeographic);
[boundaryX,boundaryY] = projfwd(targetCrs,boundaryLat,boundaryLon);
valid = isfinite(boundaryX) & isfinite(boundaryY);
if nnz(valid) < 4
    error('PHASE:InvalidProjectedExtent', ...
        'Could not project the raster extent for %s.',filePath);
end
xLimits = [min(boundaryX(valid)) max(boundaryX(valid))];
yLimits = [min(boundaryY(valid)) max(boundaryY(valid))];
if diff(xLimits) <= 0 || diff(yLimits) <= 0
    error('PHASE:InvalidProjectedExtent', ...
        'Projected raster extent is empty for %s.',filePath);
end

% Preserve approximately the original pixel count while making target cells
% close to square in the projected CRS.
rasterSize = size(source);
pixelCount = double(rasterSize(1)) * double(rasterSize(2));
aspect = diff(xLimits) / diff(yLimits);
columns = max(1,round(sqrt(pixelCount * aspect)));
rows = max(1,round(pixelCount / columns));
targetRef = maprefcells(xLimits,yLimits,[rows columns], ...
    'ColumnsStartFrom','north');

bands = max(1,size(source,3));
warped = NaN(rows,columns,bands,'single');
[targetXVector,targetYVector] = worldGrid(targetRef,'gridvectors');
blockRows = max(1,min(rows,floor(2e6 / max(1,columns))));
for firstRow = 1:blockRows:rows
    lastRow = min(rows,firstRow + blockRows - 1);
    [targetX,targetY] = meshgrid( ...
        targetXVector,targetYVector(firstRow:lastRow));
    [targetLat,targetLon] = projinv(targetCrs,targetX,targetY);
    if ~isGeographic
        [sourceX,sourceY] = projfwd(sourceCrs,targetLat,targetLon);
    end
    for band = 1:bands
        values = source(:,:,band);
        if isGeographic
            sampled = geointerp(values,sourceRef,targetLat,targetLon,'linear');
        else
            sampled = mapinterp(values,sourceRef,sourceX,sourceY,'linear');
        end
        warped(firstRow:lastRow,:,band) = single(sampled);
    end
end
if bands == 1, warped = warped(:,:,1); end
atomicGeoTiffWrite(filePath,warped,targetRef,epsgCode);
end

function [lat,lon] = rasterBoundary(reference,sourceCrs,isGeographic)
samples = linspace(0,1,101);
if isGeographic
    lonMin = reference.LongitudeLimits(1); lonMax = reference.LongitudeLimits(2);
    latMin = reference.LatitudeLimits(1); latMax = reference.LatitudeLimits(2);
    lon = [lonMin + samples*(lonMax-lonMin), ...
        repmat(lonMax,1,numel(samples)), ...
        lonMax - samples*(lonMax-lonMin), ...
        repmat(lonMin,1,numel(samples))];
    lat = [repmat(latMin,1,numel(samples)), ...
        latMin + samples*(latMax-latMin), ...
        repmat(latMax,1,numel(samples)), ...
        latMax - samples*(latMax-latMin)];
else
    xMin = reference.XWorldLimits(1); xMax = reference.XWorldLimits(2);
    yMin = reference.YWorldLimits(1); yMax = reference.YWorldLimits(2);
    x = [xMin + samples*(xMax-xMin), ...
        repmat(xMax,1,numel(samples)), ...
        xMax - samples*(xMax-xMin), ...
        repmat(xMin,1,numel(samples))];
    y = [repmat(yMin,1,numel(samples)), ...
        yMin + samples*(yMax-yMin), ...
        repmat(yMax,1,numel(samples)), ...
        yMax - samples*(yMax-yMin)];
    [lat,lon] = projinv(sourceCrs,x,y);
end
end

function atomicGeoTiffWrite(filePath,data,reference,epsgCode)
folder = fileparts(filePath);
temporary = [tempname(folder) '.tif'];
cleanup = onCleanup(@() deleteIfPresent(temporary)); %#ok<NASGU>
geotiffwrite(temporary,data,reference,'CoordRefSysCode',double(epsgCode));
[ok,message] = movefile(temporary,filePath,'f');
if ~ok
    error('PHASE:GeoTiffReplaceFailed', ...
        'Could not replace reprojected GeoTIFF %s: %s',filePath,message);
end
end

function deleteIfPresent(pathValue)
if isfile(pathValue), delete(pathValue); end
end
