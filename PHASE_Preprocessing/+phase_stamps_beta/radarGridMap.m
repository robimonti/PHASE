function frame = radarGridMap(source, request)
%RADARGRIDMAP Bounded viewport data from the same grid used by diagnostics.
% Grid curves use zero-based cell boundaries, never a shifted map overlay.
if nargin<2,request=struct();end
b=source.extent;
if isfield(request,'bounds') && ~isempty(request.bounds)
    b=double(request.bounds(:)');
end
if numel(b)~=4||any(~isfinite(b))||b(1)>=b(2)||b(3)>=b(4)|| ...
        b(1)<-180||b(2)>180||b(3)<-90||b(4)>90
    error('PHASE:GridMapBounds','Invalid geographic map bounds.');
end
stage=0;
if isfield(request,'stage'),stage=double(request.stage);end
if ~isscalar(stage)||~ismember(stage,0:3),error('PHASE:GridMapStage','Invalid candidate stage.');end
c=source.candidates;g=source.geo;
inside=isfinite(c.lon)&isfinite(c.lat)&c.lon>=b(1)&c.lon<=b(2)&c.lat>=b(3)&c.lat<=b(4);
counts=arrayfun(@(s)nnz(inside&c.stage==s),1:3);
if stage>0,inside=inside&c.stage==stage;end
ix=find(inside);total=numel(ix);limit=40000;
if total>limit,ix=ix(unique(round(linspace(1,total,limit))));end
% [longitude, latitude, D_A, stage, azimuth index, range index].
points=[c.lon(ix) c.lat(ix) c.da(ix) c.stage(ix) c.azimuth(ix) c.range(ix)];
az=[0 g.length-1];rg=[0 g.width-1];found=false;
if strcmp(g.mode,'grid')
    visible=isfinite(g.lon)&isfinite(g.lat)&g.lon>=b(1)&g.lon<=b(2)&g.lat>=b(3)&g.lat<=b(4);
    rows=find(any(visible,2));cols=find(any(visible,1));
    if ~isempty(rows)&&~isempty(cols)
        az=[rows(1)-1 rows(end)-1];rg=[cols(1)-1 cols(end)-1];found=true;
    end
else
    selected=find(inside);
    if ~isempty(selected)
        az=[min(c.azimuth(selected)) max(c.azimuth(selected))];
        rg=[min(c.range(selected)) max(c.range(selected))];found=true;
    end
end
azStep=max(1,ceil((diff(az)+1)/90));rgStep=max(1,ceil((diff(rg)+1)/90));
lines={};
if found
    azLimits=[max(-.5,az(1)-.5) min(g.length-.5,az(2)+.5)];
    rgLimits=[max(-.5,rg(1)-.5) min(g.width-.5,rg(2)+.5)];
    azEdges=unique([azLimits(1):azStep:azLimits(2) azLimits(2)]);
    rgEdges=unique([rgLimits(1):rgStep:rgLimits(2) rgLimits(2)]);
    azCurve=linspace(azLimits(1),azLimits(2),min(512,ceil(diff(azLimits))+1));
    rgCurve=linspace(rgLimits(1),rgLimits(2),min(512,ceil(diff(rgLimits))+1));
    for edge=azEdges
        lines{end+1}=curve(g,edge*ones(size(rgCurve)),rgCurve); %#ok<AGROW>
    end
    for edge=rgEdges
        lines{end+1}=curve(g,azCurve,edge*ones(size(azCurve))); %#ok<AGROW>
    end
end
frame=struct('available',true,'extent',source.extent,'bounds',b,'points',points, ...
    'lines',{lines},'stageCounts',counts,'visibleCount',total,'displayedCount',numel(ix), ...
    'azimuthStep',azStep,'rangeStep',rgStep,'runName',source.runName, ...
    'geolocationMode',source.geolocationMode,'threshold',source.threshold);
end

function xy=curve(g,az,rg)
az=az(:);rg=rg(:);
if strcmp(g.mode,'grid')
    lon=g.lonF(az+1,rg+1);lat=g.latF(az+1,rg+1);
else
    lon=g.lon(rg,az);lat=g.lat(rg,az);
end
bad=~isfinite(lon)|~isfinite(lat)|abs(lon)>180|abs(lat)>90|(lon==0&lat==0);
lon(bad)=NaN;lat(bad)=NaN;xy=[lon(:) lat(:)];
end
