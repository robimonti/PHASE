function diagnoseStampsSelectionSelfTest()
% Test spatial retention diagnostics without any real SAR/project data.
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
processing=fullfile(folder,'processing');geo=fullfile(processing,'INSAR_20200101','geo');mkdir(geo);
work=fullfile(folder,'PSI');patch=fullfile(work,'PATCH_1');mkdir(patch);
master_date='20200101';project_path=folder;psi_area_mode='full';psi_polygon=[];
save(fullfile(work,'input_StaMPS.mat'),'master_date','project_path','psi_area_mode','psi_polygon');
lon_min=9;lon_max=12;lat_min=45;lat_max=47;save(fullfile(processing,'input_preprocessing.mat'),'lon_min','lon_max','lat_min','lat_max');
lon=repmat(9:12,3,1);lat=repmat((45:47)',1,4);
lo=fullfile(geo,'20200101.lon');la=fullfile(geo,'20200101.lat');
writeRaster(lo,lon);writeRaster(la,lat);
write(fullfile(work,'psclonlat.in'),sprintf('4\n%s\n%s\n',lo,la));
write(fullfile(fileparts(geo),'phase_preprocessing_aoi.json'),jsonencode(struct('algorithm','terrain-aware-v1','requestedBounds',[9 12 45 47])));
% Exclude one valid ground pixel deliberately: the audit must find exactly one.
fid=fopen(fullfile(work,'psi_geolocation.mask'),'w');fwrite(fid,uint8([1 zeros(1,11)]),'uint8');fclose(fid);
% Initial candidates -> selection -> reselection -> retained/merged/export.
lonlat=[9 45;10 45;11 46;12 47];save(fullfile(patch,'ps1.mat'),'lonlat');
D_A=[.1;.2;.3;.4];save(fullfile(patch,'da1.mat'),'D_A');
coh_ps=[.9;.8;.7;.6];save(fullfile(patch,'pm1.mat'),'coh_ps');
ix=[1;3;4];keep_ix=[true;false;true];coh_ps2=[.95;.2;.8];coh_thresh=.6;
save(fullfile(patch,'select1.mat'),'ix','keep_ix','coh_ps2','coh_thresh');
lonlat=[12 47];save(fullfile(patch,'ps2.mat'),'lonlat');save(fullfile(work,'ps2.mat'),'lonlat');
mkdir(fullfile(work,'EXPORT'));point=1;lon=12;lat=47;writetable(table(point,lon,lat),fullfile(work,'EXPORT','displacement.csv'));
% Both SLCs finite/nonzero except one distinct missing pixel in acquisition 2.
slc=fullfile(processing,'INSAR_20200101','rslc');mkdir(slc);first=fullfile(slc,'a.rslc');second=fullfile(slc,'b.rslc');
v=ones(2,12,'single');fid=fopen(first,'w','ieee-be');fwrite(fid,v,'single');fclose(fid);
v(:,12)=0;fid=fopen(second,'w','ieee-be');fwrite(fid,v,'single');fclose(fid);
write(fullfile(work,'calamp.in'),sprintf('%s\n%s\n',first,second));
r=diagnose_stamps_grid(work);
assert(~isfield(r,'coverageAuditError'),'Coverage failed: %s',jsonencode(r));
assert(r.coverageAudit.maskAudit.unexpectedlyExcludedPixels==1);
assert(r.coverageAudit.maskAudit.unexpectedlyIncludedPixels==0);
assert(sum(r.coverageAudit.validGeolocationCounts,'all')==12);
assert(sum(r.coverageAudit.currentMaskIncludedCounts,'all')==11);
assert(sum(r.coverageAudit.commonSlcSupportCounts,'all')==11);
assert(isempty(r.selectionStageErrors));
extracted=fullfile(folder,'unpacked');unzip(fullfile(work,'PHASE_CSK_diagnostics.zip'),extracted);
a=readtable(fullfile(extracted,'selection_stages','PATCH_1_candidates.csv'));
assert(height(a)==4&&isequal(a.amplitude_dispersion,D_A));
a=readtable(fullfile(extracted,'selection_stages','PATCH_1_selected.csv'));
assert(height(a)==3&&isequal(a.longitude,[9;11;12]));
a=readtable(fullfile(extracted,'selection_stages','PATCH_1_reselected.csv'));
assert(height(a)==2&&isequal(a.longitude,[9;12]));
a=readtable(fullfile(extracted,'selection_stages','PATCH_1_retained.csv'));assert(height(a)==1&&a.longitude==12);
a=readtable(fullfile(extracted,'selection_stages','ROOT_csv_displacement.csv'));assert(height(a)==1&&a.latitude==47);
assert(isempty(dir(fullfile(extracted,'**','*.rslc'))),'SAR files must not be copied.');
% Index/mask mismatches are reported, not silently interpreted as point loss.
keep_ix=true(5,1);save(fullfile(patch,'select1.mat'),'ix','keep_ix');
r=diagnose_stamps_grid(work,fullfile(folder,'second.json'));assert(~isempty(r.selectionStageErrors));
fprintf('Selection diagnostics passed: masks, nonzero stack support, stage mapping, quality fields, CSV and safe error reporting.\n');
end
function write(file,text)
fid=fopen(file,'w');assert(fid>=0);fprintf(fid,'%s',text);fclose(fid);
end
function writeRaster(file,values)
fid=fopen(file,'w','ieee-be');fwrite(fid,values.','single');fclose(fid);
end
