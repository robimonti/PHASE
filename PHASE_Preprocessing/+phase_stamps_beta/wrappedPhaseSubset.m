function [indices, phaseRadians] = wrappedPhaseSubset(lonlat, complexPhase, center, radiusMeters)
%WRAPPEDPHASESUBSET Select raw wrapped PS phases using a geodesic radius.
% A zero radius means all persistent scatterers, matching the default UI.

if size(lonlat,2) < 2 || size(lonlat,1) ~= size(complexPhase,1)
    error('PHASE_StaMPS:wrappedShape', ...
        'Wrapped phase rows must match the StaMPS PS coordinate rows.');
end
if ~isscalar(radiusMeters) || ~isfinite(radiusMeters) || radiusMeters < 0
    error('PHASE_StaMPS:wrappedRadius', ...
        'Wrapped reference radius must be a nonnegative number of metres.');
end
if radiusMeters == 0
    indices = (1:size(lonlat,1))';
else
    center = double(center(:)');
    if numel(center) ~= 2 || any(~isfinite(center)) || ...
            abs(center(1)) > 180 || abs(center(2)) > 90
        error('PHASE_StaMPS:wrappedCenter', ...
            'Wrapped reference centre must be valid longitude and latitude.');
    end
    lat = deg2rad(double(lonlat(:,2)));
    lon = deg2rad(double(lonlat(:,1)));
    centerLat = deg2rad(center(2));
    centerLon = deg2rad(center(1));
    haversine = sin((lat-centerLat)/2).^2 + ...
        cos(lat).*cos(centerLat).*sin((lon-centerLon)/2).^2;
    distanceMeters = 2*6371008.8*asin(sqrt(min(1,max(0,haversine))));
    indices = find(distanceMeters <= radiusMeters);
end
if isempty(indices)
    error('PHASE_StaMPS:wrappedNoPoints', ...
        'No persistent scatterers fall inside the wrapped reference radius.');
end
phaseRadians = angle(complexPhase(indices,:));
end
