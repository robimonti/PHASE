function [mask,info] = preparePsiArea(cfg,workDir)
%PREPAREPSIAREA Generate StaMPS's uint8 exclusion mask in radar pixel order.
% 0 includes a pixel; 1 excludes it. Source SNAP rasters are read only.
mask='';info=struct();
isSubarea=strcmp(cfg.psi_area_mode,'subarea');
if ~isSubarea
    if isfile(fullfile(workDir,'psi_area.json'))
        error('PHASE:PsiAreaDatasetMismatch','This is a sub-area PSI dataset. Select the original full-area dataset in the hub.');
    end
    if cfg.stamps_preparation~=0
        phase_stamps_beta.validatePreparedGeolocation(workDir);
        return
    end
    polygon=zeros(0,2);prefix='psi_geolocation';
else
    polygon=phase_stamps_beta.validatePsiPolygon(cfg.psi_polygon);
    [expected,~]=phase_stamps_beta.psiAreaWorkDir(workDir,cfg);
    if ~strcmp(char(java.io.File(expected).getCanonicalPath()),char(java.io.File(workDir).getCanonicalPath()))
        error('PHASE:PsiAreaDatasetMismatch','Run the selected PSI area in its separate AOI dataset.');
    end
    prefix='psi_area';
end
grid=phase_stamps_beta.psiGrid(cfg);
mask=fullfile(workDir,[prefix '.mask']);manifest=fullfile(workDir,[prefix '.json']);
info=struct('polygon',polygon,'grid',grid.signature,'includedPixels',0,'totalPixels',grid.width*grid.height,'validPixels',0,'invalidGeolocationPixels',0);
if cfg.stamps_preparation==0 && ~strcmp(cfg.stamps_first_step,'1')
    error('PHASE:PsiAreaNeedsPreparation','Preparing new PSI candidates requires first step 1. Disable Prepare data to resume an existing area from a later step.');
end
if cfg.stamps_preparation~=0 || ~strcmp(cfg.stamps_first_step,'1')
    if ~isfile(mask)||~isfile(manifest)
        error('PHASE:PsiAreaNeedsPreparation','Prepare this PSI area from Step 1 before resuming.');
    end
    old=jsondecode(fileread(manifest));m=dir(mask);
    if ~isequal(old.polygon,polygon)||~isequaln(old.grid,grid.signature)||m.bytes~=info.totalPixels
        error('PHASE:PsiAreaChanged','PSI area or SNAP geolocation changed. Enable Prepare data and restart at Step 1.');
    end
    info=old;return
end
if isSubarea
    p=polygon;p(:,1)=rad2deg(unwrap(deg2rad(p(:,1))));centre=mean(p(1:end-1,1));
end
temporary=[tempname(workDir) '.mask'];
out=fopen(temporary,'w');lo=fopen(grid.longitude,'r',grid.byteOrder);la=fopen(grid.latitude,'r',grid.byteOrder);
cleanup=onCleanup(@()closeFiles([out lo la],temporary));
if any([out lo la]<0),error('PHASE:PsiMaskIO','Cannot open PSI mask or SNAP geolocation files.');end
remaining=info.totalPixels;
while remaining>0
    count=min(remaining,1000000);
    lon=fread(lo,count,'single=>double');lat=fread(la,count,'single=>double');
    if numel(lon)~=count||numel(lat)~=count,error('PHASE:PsiGridReadFailed','Incomplete SNAP geolocation raster.');end
    valid=isfinite(lon)&isfinite(lat)&abs(lon)<=180&abs(lat)<90&~(lon==0&lat==0);
    info.validPixels=info.validPixels+nnz(valid);
    inside=valid;
    if isSubarea
        lon=centre+mod(lon-centre+180,360)-180;
        inside=inside&inpolygon(lon,lat,p(:,1),p(:,2));
    end
    if fwrite(out,uint8(~inside),'uint8')~=count,error('PHASE:PsiMaskIO','Incomplete mask write. Check free disk space.');end
    info.includedPixels=info.includedPixels+nnz(inside);remaining=remaining-count;
end
info.invalidGeolocationPixels=info.totalPixels-info.validPixels;
if info.includedPixels==0
    if ~isSubarea,error('PHASE:NoValidGeolocation','The SNAP export has no valid radar geolocation. Check the DEM and preprocessing geolocation outputs.');end
    error('PHASE:PsiAreaEmpty','The polygon contains no valid SNAP radar pixels. Draw an area inside the dataset coverage.');
end
fclose(out);
movefile(temporary,mask,'f');
clear cleanup % Close readers before any other fopen can reuse their IDs.
fid=fopen(manifest,'w');if fid<0,error('PHASE:PsiMaskIO','Cannot save PSI area metadata.');end
fwrite(fid,jsonencode(info),'char');fclose(fid);
if ~isSubarea && info.invalidGeolocationPixels==0
    delete(mask);mask=''; % Preserve the original command for fully valid stacks.
end
end

function closeFiles(handles,temporary)
for handle=handles
    if handle>0,try,fclose(handle);catch,end,end
end
if isfile(temporary),delete(temporary);end
end
