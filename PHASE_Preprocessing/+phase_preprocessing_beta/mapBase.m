function lines = mapBase()
%MAPBASE Local coastline vectors for the dependency-free HTML map.

persistent cached
if ~isempty(cached), lines = cached; return; end
lines = repmat(struct('coordinates', zeros(0, 2)), 0, 1);
try
    data = load('coastlines');
    lat = double(data.coastlat(:)); lon = double(data.coastlon(:));
    boundaries = [0; find(isnan(lat) | isnan(lon)); numel(lat) + 1];
    for k = 1:numel(boundaries)-1
        first = boundaries(k) + 1; last = boundaries(k+1) - 1;
        if last - first < 1, continue; end
        coords = [lon(first:last), lat(first:last)];
        if size(coords, 1) > 500
            stride = ceil(size(coords, 1) / 500);
            coords = coords([1:stride:end, end], :);
        end
        lines(end+1).coordinates = coords; %#ok<AGROW>
    end
catch
    % The map still provides the graticule, AOI and footprints if the
    % coastline dataset is unavailable.
end
cached = lines;
end
