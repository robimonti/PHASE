function grid = psiGrid(cfg)
%PSIGRID Match mt_prep_snap's master parameter and geolocation selection.
root=phase_stamps_beta.exportPath(cfg.project_path,cfg.master_date);
pars=dir(fullfile(root,'*slc',[cfg.master_date '.*slc.par']));
lons=dir(fullfile(root,'geo','*.lon'));lats=dir(fullfile(root,'geo','*.lat'));
if isempty(pars)||isempty(lons)||isempty(lats)
    error('PHASE:PsiGridMissing','SNAP master parameters and geo/*.lon/*.lat are required in %s.',root);
end
pars=sortEntries(pars);lons=sortEntries(lons);lats=sortEntries(lats);
parPath=fullfile(pars(1).folder,pars(1).name);text=fileread(parPath);
w=regexp(text,'(?m)^\s*range_samples\s*:?\s*(\d+)','tokens','once');
h=regexp(text,'(?m)^\s*azimuth_lines\s*:?\s*(\d+)','tokens','once');
if isempty(w)||isempty(h),error('PHASE:PsiGridInvalid','Missing radar dimensions in %s.',parPath);end
grid=struct('width',str2double(w{1}),'height',str2double(h{1}), ...
    'longitude',fullfile(lons(1).folder,lons(1).name), ...
    'latitude',fullfile(lats(1).folder,lats(1).name),'source',root);
if grid.width<1||grid.height<1||lons(1).bytes~=4*grid.width*grid.height||lats(1).bytes~=lons(1).bytes
    error('PHASE:PsiGridInvalid','SNAP geolocation sizes do not match the master radar dimensions.');
end
% SNAP/StaMPS float32 rasters are big endian (mt_prep_snap uses byteswap=1
% on the supported little-endian platforms). Do not guess byte order.
grid.byteOrder='ieee-be';
grid.signature=struct('source',char(java.io.File(root).getCanonicalPath()),'master',cfg.master_date,'width',grid.width,'height',grid.height, ...
    'longitude',lons(1).name,'latitude',lats(1).name, ...
    'lonBytes',lons(1).bytes,'latBytes',lats(1).bytes, ...
    'lonModified',lons(1).datenum,'latModified',lats(1).datenum, ...
    'parModified',pars(1).datenum);
[~,area]=phase_stamps_beta.preprocessingArea(root);
if ~isempty(area),grid.signature.preprocessingAOI=area;end
end

function out=sortEntries(out)
paths=arrayfun(@(entry)fullfile(entry.folder,entry.name),out,'UniformOutput',false);
[~,order]=sort(paths);out=out(order);
end
