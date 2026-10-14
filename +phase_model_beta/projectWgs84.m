function [x,y] = projectWgs84(crs,lonlat)
%PROJECTWGS84 Project geographic coordinates in the selected project CRS.
if size(lonlat,2)~=2 || any(~isfinite(lonlat(:)))
    error('PHASE_Model_beta:invalidCoordinates', ...
        'Project coordinates must be finite longitude/latitude pairs.');
end
[x,y] = projfwd(crs,lonlat(:,2),lonlat(:,1));
if any(~isfinite(x)) || any(~isfinite(y))
    error('PHASE_Model_beta:projectionFailed', ...
        'The selected CRS could not project all coordinates.');
end
end
