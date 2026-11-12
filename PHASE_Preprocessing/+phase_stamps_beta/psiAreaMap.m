function data = psiAreaMap(cfg)
%PSIAREAMAP Lightweight radar-coverage preview, never used as the mask.
grid=phase_stamps_beta.psiGrid(cfg);
rows=unique(round(linspace(1,grid.height,min(64,grid.height))));
cols=unique(round(linspace(1,grid.width,min(64,grid.width))));
lon=sample(grid.longitude);lat=sample(grid.latitude);
valid=isfinite(lon)&isfinite(lat)&abs(lon)<=180&abs(lat)<90&~(lon==0&lat==0);
[bounds,~]=phase_stamps_beta.preprocessingArea(grid.source);
if ~isempty(bounds),valid=valid&lon>=bounds(1)&lon<=bounds(2)&lat>=bounds(3)&lat<=bounds(4);end
points=unique([lon(valid),lat(valid)],'rows');
if size(points,1)<3,error('PHASE:PsiGridInvalid','No valid radar coverage found in SNAP geolocation rasters.');end
points(:,1)=rad2deg(unwrap(deg2rad(points(:,1))));
ring=points(convhull(points(:,1),points(:,2)),:);
ring(:,1)=mod(ring(:,1)+180,360)-180;
footprint=struct('name','SNAP radar coverage (preview)','coordinates',ring, ...
    'source','SNAP export','selected',true);
data=struct('footprints',{{footprint}},'coastlines',{{}}, ...
    'width',grid.width,'height',grid.height,'message', ...
    'The outline previews valid radar coverage. Only valid pixels inside your polygon are selected.');
if ~isempty(bounds)
    ring=[bounds(1) bounds(3);bounds(2) bounds(3);bounds(2) bounds(4);bounds(1) bounds(4);bounds(1) bounds(3)];
    data.footprints{end+1}=struct('name','Requested preprocessing ground AOI','coordinates',ring,'source','PHASE preprocessing','selected',false);
    data.message='Coverage is shown inside the requested ground AOI. The registration halo is excluded from PSI.';
end

    function values=sample(path)
        fid=fopen(path,'r',grid.byteOrder);
        if fid<0,error('PHASE:PsiGridReadFailed','Cannot read %s.',path);end
        cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
        values=zeros(numel(rows),numel(cols));
        for k=1:numel(rows)
            fseek(fid,4*(rows(k)-1)*grid.width,'bof');
            line=fread(fid,grid.width,'single=>double');values(k,:)=line(cols);
        end
        values=values(:);
    end
end
