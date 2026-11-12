function psiAreaSelfTest()
%PSIAREASELFTEST Radar mask alignment, dataset isolation and resume safeguards.
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
cfg=phase_stamps_beta.defaultConfig();cfg.project_path=folder;cfg.master_date='20200101';
root=phase_stamps_beta.exportPath(folder,cfg.master_date);
mkdir(fullfile(root,'rslc'));mkdir(fullfile(root,'geo'));
writeText(fullfile(root,'rslc','20200101.rslc.par'),'range_samples: 4\nazimuth_lines: 3\n');
lon=[9 10 11 12;9 10 11 12;9 10 11 12];lat=[45 45 45 45;46 46 46 46;47 47 47 47];
writeRaster('lon',lon);writeRaster('lat',lat);
base=fullfile(folder,'ASC_test');mkdir(base);
writeText(fullfile(base,'original.txt'),'Keep the full dataset.');
assert(isempty(phase_stamps_beta.preparePsiArea(cfg,base)));
% Full-area processing must mask invalid coordinates, too (not just polygons).
badLon=lon;badLat=lat;badLon(1,1)=0;badLat(1,1)=0;badLon(1,2)=NaN;
writeRaster('lon',badLon);writeRaster('lat',badLat);
[fullMask,fullInfo]=phase_stamps_beta.preparePsiArea(cfg,base);
assert(isfile(fullMask)&&fullInfo.includedPixels==10&&fullInfo.invalidGeolocationPixels==2);
fid=fopen(fullMask,'r');bytes=fread(fid,Inf,'uint8');fclose(fid);
assert(isequal(bytes(:).',[1 1 zeros(1,10)]));
lonlat=[0 0;9 45];save(fullfile(base,'ps1.mat'),'lonlat');
resumeCfg=cfg;resumeCfg.stamps_preparation=1;
mustFail(@()phase_stamps_beta.preparePsiArea(resumeCfg,base),'PHASE:InvalidPreparedGeolocation');
delete(fullfile(base,'ps1.mat'));
writeRaster('lon',lon);writeRaster('lat',lat);
cfg.psi_area_mode='subarea';cfg.psi_polygon=[9.5 45.5;11.5 45.5;11.5 47.5;9.5 47.5;9.5 45.5];
[target,cfg]=phase_stamps_beta.psiAreaWorkDir(base,cfg);mkdir(target);
assert(~strcmp(base,target)&&strcmp(cfg.psi_source_dataset,'ASC_test'));
phase_stamps_beta.saveConfig(target,cfg);
[mask,info]=phase_stamps_beta.preparePsiArea(cfg,target);
fid=fopen(mask,'r');actual=fread(fid,Inf,'uint8=>uint8');fclose(fid);
expected=uint8([1 1 1 1;1 0 0 1;1 0 0 1]);
assert(isequal(actual,reshape(expected.',[],1)),'Radar row/column mask alignment changed.');
assert(info.includedPixels==4&&info.totalPixels==12);
assert(isequal(readRaster('lon'),lon)&&isfile(fullfile(base,'original.txt')));
preview=phase_stamps_beta.psiAreaMap(cfg);assert(numel(preview.footprints)==1);
assert(size(preview.footprints{1}.coordinates,2)==2);
ui=phase_stamps_beta.configToUi(cfg);roundtrip=phase_stamps_beta.uiToConfig(ui,cfg);
assert(phase_stamps_beta.configsEqual(cfg,roundtrip));
loaded=phase_stamps_beta.loadConfig(target);assert(isequal(cfg.psi_polygon,loaded.psi_polygon));
cfg.stamps_preparation=1;mkdir(fullfile(target,'PATCH_1'));
[resumed,resumedInfo]=phase_stamps_beta.preparePsiArea(cfg,target);
assert(strcmp(mask,resumed)&&resumedInfo.includedPixels==4);
unsafe=cfg;unsafe.stamps_preparation=0;unsafe.stamps_first_step='8';
mustFail(@()phase_stamps_beta.preparePsiArea(unsafe,target),'PHASE:PsiAreaNeedsPreparation');
changed=cfg;changed.psi_polygon(:,1)=changed.psi_polygon(:,1)+0.1;
mustFail(@()phase_stamps_beta.psiAreaWorkDir(target,changed),'PHASE:PsiAreaNeedsPreparation');
cfg.psi_area_mode='full';[full,cfgFull]=phase_stamps_beta.psiAreaWorkDir(target,cfg);
assert(strcmp(full,base)&&isempty(cfgFull.psi_source_dataset));
cfg.psi_area_mode='subarea';
% Same dimensions but changed geolocation must never silently reuse candidates.
writeRaster('lon',lon+0.2);
changedFile=java.io.File(fullfile(root,'geo','20200101.lon'));
assert(changedFile.setLastModified(changedFile.lastModified()+2000));
mustFail(@()phase_stamps_beta.preparePsiArea(cfg,target),'PHASE:PsiAreaChanged');
writeRaster('lon',lon);cfg.stamps_preparation=0;
% Invalid geolocation pixels are always excluded, including polygon boundaries.
bad=lon;bad(2,2)=NaN;writeRaster('lon',bad);
[~,info]=phase_stamps_beta.preparePsiArea(cfg,target);assert(info.includedPixels==3);
cfg.psi_polygon=[70 10;71 10;71 11;70 11;70 10];
[emptyTarget,cfg]=phase_stamps_beta.psiAreaWorkDir(target,cfg);mkdir(emptyTarget);
mustFail(@()phase_stamps_beta.preparePsiArea(cfg,emptyTarget),'PHASE:PsiAreaEmpty');
mustFail(@()phase_stamps_beta.validatePsiPolygon([0 0;1 1;0 1;1 0;0 0]),'PHASE:InvalidPsiArea');
% Antimeridian selection: +179 and -179 must be neighbours in radar space.
writeRaster('lon',[178 179 -179 -178;178 179 -179 -178;178 179 -179 -178]);
cfg.psi_polygon=[178.5 44.5;-178.5 44.5;-178.5 47.5;178.5 47.5;178.5 44.5];
[dateTarget,cfg]=phase_stamps_beta.psiAreaWorkDir(target,cfg);mkdir(dateTarget);
[~,info]=phase_stamps_beta.preparePsiArea(cfg,dateTarget);assert(info.includedPixels==6);
% Terrain-aware CSK halos are excluded even in full-area PSI mode.
writeRaster('lon',lon);writeRaster('lat',lat);
area=struct('algorithm','terrain-aware-v1','requestedBounds',[9.5 11.5 45.5 47.5]);
areaFile=fullfile(root,'phase_preprocessing_aoi.json');writeText(areaFile,jsonencode(area));
cfg.psi_area_mode='full';cfg.psi_source_dataset='';cfg.stamps_preparation=0;
[groundMask,groundInfo]=phase_stamps_beta.preparePsiArea(cfg,base);
assert(isfile(groundMask)&&groundInfo.includedPixels==4&&groundInfo.outsidePreprocessingAOIPixels==8);
assert(groundInfo.invalidGeolocationPixels==0);
preview=phase_stamps_beta.psiAreaMap(cfg);assert(numel(preview.footprints)==2);
cfg.stamps_preparation=1;phase_stamps_beta.preparePsiArea(cfg,base);
area.requestedBounds=[10.5 11.5 45.5 47.5];writeText(areaFile,jsonencode(area));
mustFail(@()phase_stamps_beta.preparePsiArea(cfg,base),'PHASE:PsiAreaChanged');
delete(areaFile);constellation='CSK';save(fullfile(fileparts(root),'input_preprocessing.mat'),'constellation');
mustFail(@()phase_stamps_beta.psiAreaMap(cfg),'PHASE:CskTerrainRebuildRequired');
fprintf('PSI area regression passed: byte order, masks, ground AOI/halo, legacy guard, isolation, resume, invalid pixels and date line.\n');

    function writeRaster(extension,values)
        fid=fopen(fullfile(root,'geo',['20200101.' extension]),'w','ieee-be');
        fwrite(fid,values.','single');fclose(fid);
    end
    function values=readRaster(extension)
        fid=fopen(fullfile(root,'geo',['20200101.' extension]),'r','ieee-be');
        values=fread(fid,[4 3],'single=>double').';fclose(fid);
    end
end

function mustFail(action,id)
try,action();catch ME,assert(strcmp(ME.identifier,id),'Unexpected error: %s',ME.message);return;end
error('PHASE:TestExpectedFailure','Expected %s.',id);
end

function writeText(path,text)
fid=fopen(path,'w');fprintf(fid,text);fclose(fid);
end
