function [crs,utmZone,epsgCode] = selectProjectUtm(lonlatAOI)
%SELECTPROJECTUTM Choose one metric CRS for the complete modelling domain.
% Adjacent UTM zones are projected into the zone covering most AOI area.

valid = all(isfinite(lonlatAOI),2);
vertices = double(lonlatAOI(valid,:));
if size(vertices,1)<3 || any(abs(vertices(:,1))>180) || ...
        any(abs(vertices(:,2))>90)
    error('PHASE_Model_beta:invalidAoiCoordinates', ...
        'The AOI must have at least three finite WGS84 longitude/latitude vertices.');
end
if range(vertices(:,1))>180
    error('PHASE_Model_beta:antimeridianAoi', ...
        'AOIs crossing the antimeridian are not supported by the current UTM workflow.');
end
if min(vertices(:,2))<0 && max(vertices(:,2))>=0
    error('PHASE_Model_beta:equatorAoi', ...
        'AOIs crossing the equator need a different projected CRS.');
end
if min(vertices(:,2))< -80 || max(vertices(:,2))>84
    error('PHASE_Model_beta:polarAoi', ...
        'UTM is not defined for this polar AOI.');
end

% Approximate the area in each zone by regular sampling within the AOI.
% cos(latitude) compensates for the varying area of longitude cells.
longitude = linspace(min(vertices(:,1)),max(vertices(:,1)),101);
latitude = linspace(min(vertices(:,2)),max(vertices(:,2)),101);
[sampleLon,sampleLat] = meshgrid(longitude,latitude);
inside = inpolygon(sampleLon,sampleLat,lonlatAOI(:,1),lonlatAOI(:,2));
if ~any(inside(:))
    error('PHASE_Model_beta:emptyAoiArea', ...
        'The AOI has no measurable interior. Check its polygon geometry.');
end
zones = min(60,max(1,floor((sampleLon(inside)+180)/6)+1));
weights = cosd(sampleLat(inside));
zoneArea = accumarray(zones(:),weights(:),[60 1],@sum,0);
[~,zone] = max(zoneArea);

vertexZones = min(60,max(1,floor((vertices(:,1)+180)/6)+1));
if any(abs(vertexZones-zone)>1)
    error('PHASE_Model_beta:wideAoi', ...
        ['The AOI extends beyond the UTM zones adjacent to the dominant zone. ', ...
         'A local projected CRS is required for this extent.']);
end

hemisphereNorth = mean(vertices(:,2))>=0;
epsgCode = (32600*hemisphereNorth + 32700*~hemisphereNorth) + zone;
crs = projcrs(epsgCode);
bands = 'CDEFGHJKLMNPQRSTUVWX';
bandIndex = min(numel(bands),max(1,floor((mean(vertices(:,2))+80)/8)+1));
utmZone = sprintf('%02d %c',zone,bands(bandIndex));
end
