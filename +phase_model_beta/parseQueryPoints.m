function points=parseQueryPoints(raw)
%PARSEQUERYPOINTS Validate manually entered station IDs and WGS84 coordinates.
try
    entries=jsondecode(char(string(raw)));
catch
    error('PHASE_Model_beta:queryPoints','Query points must be a valid list of ID, longitude and latitude rows.');
end
if isempty(entries)||~isstruct(entries)||~all(isfield(entries,{'id','longitude','latitude'}))
    error('PHASE_Model_beta:queryPoints','Add at least one query point with ID, longitude and latitude.');
end
ids=strings(numel(entries),1);lon=zeros(numel(entries),1);lat=lon;
for k=1:numel(entries)
    id=entries(k).id;
    if ~(ischar(id)||(isstring(id)&&isscalar(id)))
        error('PHASE_Model_beta:queryPoints','Query-point IDs must be text.');
    end
    ids(k)=strtrim(string(id));
    if strlength(ids(k))==0 || ~isempty(regexp(char(ids(k)),'[<>:"/\\|?*\x00-\x1f]','once')) || ...
            endsWith(ids(k),'.') || endsWith(ids(k),' ')
        error('PHASE_Model_beta:queryPoints','Use a non-empty point ID without filename separators or reserved characters.');
    end
    lonValue=entries(k).longitude;latValue=entries(k).latitude;
    if ~isnumeric(lonValue)||~isscalar(lonValue)||~isfinite(lonValue)||abs(lonValue)>180 || ...
            ~isnumeric(latValue)||~isscalar(latValue)||~isfinite(latValue)||abs(latValue)>90
        error('PHASE_Model_beta:queryPoints','Point %s needs finite WGS84 longitude (-180 to 180) and latitude (-90 to 90).',ids(k));
    end
    lon(k)=lonValue;lat(k)=latValue;
end
if numel(unique(lower(ids)))~=numel(ids)
    error('PHASE_Model_beta:queryPoints','Query-point IDs must be unique (case insensitive).');
end
points=table(ids,lon,lat,'VariableNames',{'id','longitude','latitude'});
end
