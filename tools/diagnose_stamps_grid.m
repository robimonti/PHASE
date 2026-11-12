function report = diagnose_stamps_grid(workDir,outputPath)
%DIAGNOSE_STAMPS_GRID Collect small, read-only PSI/geolocation diagnostics.
% Select a PSI dataset containing input_StaMPS.mat and/or PATCH_*.
% Only the JSON report is written; no inputs, parameters or results change.
if nargin<1 || isempty(workDir)
    workDir=uigetdir(pwd,'Select the PSI dataset folder (containing PATCH_*)');
    if isequal(workDir,0),report=struct();return;end
end
workDir=char(java.io.File(workDir).getCanonicalPath());
selectedDir=workDir;
[~,selectedName]=fileparts(workDir);
if startsWith(upper(selectedName),'PATCH_') && isfile(fullfile(fileparts(workDir),'input_StaMPS.mat'))
    workDir=fileparts(workDir);
end
if ~isfile(fullfile(workDir,'input_StaMPS.mat'))
    error('PHASE:DiagnosticDataset','Select the PSI dataset containing input_StaMPS.mat, or its PATCH_1 folder.');
end
if nargin<2,outputPath=fullfile(selectedDir,'PHASE_PSI_grid_diagnostics.json');end
report=struct('diagnosticSchema',2,'generatedAt',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'dataset',workDir,'config',struct(),'patches',{{}},'geolocation',struct());
configPath=fullfile(workDir,'input_StaMPS.mat');
if isfile(configPath)
    names={'master_date','filter_grid_size','amplitude_threshold','stamps_preparation', ...
        'stamps_first_step','stamps_last_step','psi_area_mode','psi_polygon', ...
        'project_path','merge_resample_size','merge_standard_dev','weed_standard_dev', ...
        'weed_neighbours','weed_zero_elevation','select_method','percent_rand','density_rand'};
    available=whos('-file',configPath);names=intersect(names,{available.name},'stable');
    report.config=load(configPath,names{:});
end
patches=dir(fullfile(workDir,'PATCH_*'));patches=patches([patches.isdir]);
for k=1:numel(patches)
    folder=fullfile(patches(k).folder,patches(k).name);
    item=struct('name',patches(k).name,'files',struct(),'savedArrays',{{}});
    psFile=fullfile(folder,'ps1.mat');
    if isfile(psFile)
        vars=whos('-file',psFile);
        item.savedArrays=arrayfun(@(v)struct('name',v.name,'size',v.size,'class',v.class), ...
            vars,'UniformOutput',false);
        requested=intersect({'lonlat','xy','n_ps','n_ifg','n_image'}, {vars.name},'stable');
        ps=load(psFile,requested{:});
        for field={'n_ps','n_ifg','n_image'}
            if isfield(ps,field{1}),item.(field{1})=ps.(field{1});end
        end
        if isfield(ps,'lonlat'),item.lonlat=coordinateStats(ps.lonlat);end
        if isfield(ps,'xy') && size(ps.xy,2)>=3
            xy=double(ps.xy(:,2:3));valid=all(isfinite(xy),2);
            item.xy=coordinateStats(xy,false);
            if any(valid) && isfield(report.config,'filter_grid_size')
                step=report.config.filter_grid_size;
                if isscalar(step)&&isfinite(step)&&step>0
                    span=max(xy(valid,:),[],1)-min(xy(valid,:),[],1);
                    cells=max(1,ceil((span+1e-6)/step)-1);
                    item.estimatedGammaGridColumnsRows=cells;
                    if isfield(ps,'n_ifg')
                        item.oneGammaArrayGiB=prod(cells)*double(ps.n_ifg)*4/1024^3;
                    end
                end
            end
        end
    end
    for suffix={'ij','ll','ph','da','hgt'}
        path=fullfile(folder,['pscands.1.' suffix{1}]);info=dir(path);
        if ~isempty(info),item.files.(suffix{1})=info.bytes;end
    end
    ll=fullfile(folder,'pscands.1.ll');
    if isfile(ll)
        item.rawLonlatBigEndian=coordinateStats(samplePairs(ll,'ieee-be'));
        item.rawLonlatLittleEndian=coordinateStats(samplePairs(ll,'ieee-le'));
    end
    item.stages=stageCounts(folder);
    report.patches{end+1}=item;
end
report.rootStages=stageCounts(workDir);
list=fullfile(workDir,'patch.list');
if isfile(list)
    lines=strtrim(splitlines(string(fileread(list))));lines(lines=="")=[];
    report.mergePatchList=cellstr(lines);
end
geoInput=fullfile(workDir,'psclonlat.in');
if isfile(geoInput)
    lines=splitlines(string(fileread(geoInput)));lines=strtrim(lines);lines(lines=="")=[];
    if numel(lines)>=3
        lonPath=strip(char(lines(2)),'both','"');latPath=strip(char(lines(3)),'both','"');
        if isfile(lonPath)&&isfile(latPath)
            report.geolocation.width=str2double(lines(1));
            report.geolocation.longitudeBytes=dirBytes(lonPath);
            report.geolocation.latitudeBytes=dirBytes(latPath);
            report.geolocation.sampleBigEndian=coordinateStats(sampleGeo(lonPath,latPath,'ieee-be'));
            report.geolocation.sampleLittleEndian=coordinateStats(sampleGeo(lonPath,latPath,'ieee-le'));
        end
    end
end
area=fullfile(workDir,'psi_area.json');
if isfile(area)
    try,report.area=jsondecode(fileread(area));catch ME,report.areaReadError=ME.message;end
end
try
    report.processingAudit=processingAudit(workDir,report);
catch ME
    report.processingAuditError=ME.message;
end
fid=fopen(outputPath,'w','n','UTF-8');
if fid<0,error('PHASE:DiagnosticWriteFailed','Cannot write %s.',outputPath);end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
fprintf('Diagnostic report saved: %s\nNo processing files or parameters were modified.\n',outputPath);
clear cleanup
if nargin<2
    bundle=writeDiagnosticBundle(workDir,report,outputPath);
    fprintf('Send this single ZIP: %s\n',bundle);
    if ispc,winopen(fileparts(bundle));end
end
end

function out=coordinateStats(values,checkBounds)
if nargin<2,checkBounds=true;end
values=double(values);valid=all(isfinite(values),2);
out=struct('rows',size(values,1),'nonfiniteRows',nnz(~valid));
if ~isempty(values)
    out.zeroPairs=nnz(all(values==0,2));
    if checkBounds,out.outsideLonLatBounds=nnz(valid & (abs(values(:,1))>180 | abs(values(:,2))>90));end
    if any(valid)
        out.minimum=min(values(valid,:),[],1);out.maximum=max(values(valid,:),[],1);
        out.median=median(values(valid,:),1);
    end
end
end

function values=samplePairs(path,order)
count=floor(dirBytes(path)/8);indices=unique(round(linspace(0,max(0,count-1),min(5000,count))));
fid=fopen(path,'r',order);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
values=nan(numel(indices),2);
for k=1:numel(indices)
    fseek(fid,8*indices(k),'bof');pair=fread(fid,2,'single=>double');
    if numel(pair)==2,values(k,:)=pair.';end
end
end

function values=sampleGeo(lonPath,latPath,order)
count=floor(min(dirBytes(lonPath),dirBytes(latPath))/4);
indices=unique(round(linspace(0,max(0,count-1),min(5000,count))));
lo=fopen(lonPath,'r',order);la=fopen(latPath,'r',order);
cleanup=onCleanup(@()closePair(lo,la)); %#ok<NASGU>
values=nan(numel(indices),2);
for k=1:numel(indices)
    fseek(lo,4*indices(k),'bof');fseek(la,4*indices(k),'bof');
    x=fread(lo,1,'single=>double');y=fread(la,1,'single=>double');
    if ~isempty(x)&&~isempty(y),values(k,:)=[x y];end
end
end

function value=dirBytes(path)
entry=dir(path);value=entry(1).bytes;
end

function closePair(a,b)
fclose(a);fclose(b);
end

function out=stageCounts(folder)
% Load only small selection masks/counts, never the full phase matrices.
out=struct();
for name={'ps1','ps2'}
    file=fullfile(folder,[name{1} '.mat']);
    if ~isfile(file),continue;end
    vars=whos('-file',file);fields=intersect({'n_ps','lonlat'},{vars.name},'stable');
    data=load(file,fields{:});item=struct();
    if isfield(data,'n_ps'),item.n_ps=data.n_ps;end
    if isfield(data,'lonlat'),item.lonlat=coordinateStats(data.lonlat);end
    out.(name{1})=item;
end
file=fullfile(folder,'select1.mat');
if isfile(file)
    vars=whos('-file',file);fields=intersect({'ix','keep_ix'},{vars.name},'stable');
    data=load(file,fields{:});item=struct();
    if isfield(data,'ix'),item.selectedBeforeReselection=numel(data.ix);end
    if isfield(data,'keep_ix'),item.keptAfterReselection=nnz(data.keep_ix);end
    out.selection=item;
end
file=fullfile(folder,'weed1.mat');
if isfile(file)
    vars=whos('-file',file);fields=intersect({'ix_weed','ix_weed2'},{vars.name},'stable');
    data=load(file,fields{:});item=struct();
    for k=1:numel(fields),item.(fields{k})=nnz(data.(fields{k}));end
    out.weeding=item;
end
file=fullfile(folder,'no_ps_info.mat');
if isfile(file),out.noPS=load(file,'stamps_step_no_ps');end
end

function out=processingAudit(workDir,report)
out=struct();
% Resolve the actual export from the paths consumed by StaMPS, rather than
% assuming a project layout or requiring the user to select another folder.
geoInput=fullfile(workDir,'psclonlat.in');
if ~isfile(geoInput),out.message='No psclonlat.in found in this PSI dataset.';return;end
lines=strtrim(splitlines(string(fileread(geoInput))));lines(lines=="")=[];
if numel(lines)<3,out.message='Incomplete psclonlat.in.';return;end
lonPath=strip(char(lines(2)),'both','"');latPath=strip(char(lines(3)),'both','"');
preprocessing=fileparts(fileparts(fileparts(lonPath)));
out.preprocessingFolder=preprocessing;
config=fullfile(preprocessing,'input_preprocessing.mat');cfg=struct();
if isfile(config)
    vars=whos('-file',config);
    names=intersect({'constellation','first_step','master_processing','master_date', ...
        'lon_min','lon_max','lat_min','lat_max','dem_name','dem_file','num_gcp'}, ...
        {vars.name},'stable');
    cfg=load(config,names{:});out.preprocessingConfig=cfg;
end
for name={'CSK_subsetMgraph2run.xml','CSK_subsetgraph2run.xml','CSK_coreg_ifg2run.xml','CSK_export2run.xml'}
    file=fullfile(preprocessing,'graphs',name{1});
    if isfile(file)
        content=fileread(file);key=matlab.lang.makeValidName(name{1});
        out.generatedGraphs.(key)=struct('modified',dirDate(file),'xml',content);
    end
end
for name={'project_master.conf','project.conf'}
    file=fullfile(preprocessing,'snap2stamps','bin',name{1});
    if isfile(file)
        key=matlab.lang.makeValidName(name{1});
        out.generatedConfigs.(key)=struct('modified',dirDate(file),'text',fileread(file));
    end
end
if isfield(report.config,'master_date')
    file=fullfile(preprocessing,'master',[char(string(report.config.master_date)) '_sub.dim']);
    if isfile(file),out.masterHeader=struct('modified',dirDate(file),'xml',fileread(file));end
end
if isfile(lonPath)&&isfile(latPath)
    out.fullExportGeolocation=scanGeo(lonPath,latPath,cfg);
end
if isfield(cfg,'dem_file')&&isfile(cfg.dem_file)
    try
        info=georasterinfo(cfg.dem_file);
        dem=struct('rasterSize',info.RasterSize,'missingDataIndicator',info.MissingDataIndicator);
        ref=info.RasterReference;
        if ~isempty(ref)
            for field={'LatitudeLimits','LongitudeLimits','XWorldLimits','YWorldLimits','RasterInterpretation'}
                if isprop(ref,field{1}),dem.(field{1})=ref.(field{1});end
            end
        end
        if ~isempty(info.CoordinateReferenceSystem)
            dem.crsName=info.CoordinateReferenceSystem.Name;
            try,dem.crsWKT=wktstring(info.CoordinateReferenceSystem);catch,end
        end
        out.dem=dem;
    catch ME,out.demReadError=ME.message;end
end
end

function out=scanGeo(lonPath,latPath,cfg)
% Full streaming scan, including nodata counts; no SAR phase data is read.
lo=fopen(lonPath,'r','ieee-be');la=fopen(latPath,'r','ieee-be');
cleanup=onCleanup(@()closePair(lo,la)); %#ok<NASGU>
if lo<0||la<0,error('PHASE:DiagnosticGeoRead','Cannot read geolocation rasters.');end
out=struct('totalPixels',0,'validPixels',0,'zeroPairs',0,'nonfinitePixels',0, ...
    'outOfBoundsPixels',0,'validInsideConfiguredAOI',0,'minimum',[Inf Inf],'maximum',[-Inf -Inf]);
box=all(isfield(cfg,{'lon_min','lon_max','lat_min','lat_max'}));
while true
    x=fread(lo,1000000,'single=>double');y=fread(la,1000000,'single=>double');
    if isempty(x)&&isempty(y),break;end
    if numel(x)~=numel(y),error('PHASE:DiagnosticGeoSize','Geolocation rasters have different lengths.');end
    finite=isfinite(x)&isfinite(y);bounds=abs(x)<=180&abs(y)<=90;zero=x==0&y==0;
    valid=finite&bounds&~zero;
    out.totalPixels=out.totalPixels+numel(x);out.validPixels=out.validPixels+nnz(valid);
    out.zeroPairs=out.zeroPairs+nnz(zero);out.nonfinitePixels=out.nonfinitePixels+nnz(~finite);
    out.outOfBoundsPixels=out.outOfBoundsPixels+nnz(finite&~bounds);
    if any(valid)
        out.minimum=min(out.minimum,[min(x(valid)),min(y(valid))]);
        out.maximum=max(out.maximum,[max(x(valid)),max(y(valid))]);
    end
    if box
        inside=valid&x>=cfg.lon_min&x<=cfg.lon_max&y>=cfg.lat_min&y<=cfg.lat_max;
        out.validInsideConfiguredAOI=out.validInsideConfiguredAOI+nnz(inside);
    end
end
out.validOutsideConfiguredAOI=out.validPixels-out.validInsideConfiguredAOI;
end

function value=dirDate(file)
entry=dir(file);value=entry(1).date;
end

function destination=writeDiagnosticBundle(workDir,report,reportPath)
% Package only geometry, metadata and logs; no SAR amplitudes or phases.
staging=tempname;mkdir(staging);cleanup=onCleanup(@()rmdir(staging,'s')); %#ok<NASGU>
copy(reportPath,'PHASE_PSI_grid_diagnostics.json');
for name={'psi_area.json','psi_geolocation.json','psclonlat.in'}
    copy(fullfile(workDir,name{1}),fullfile('PSI',name{1}));
end
geoInput=fullfile(workDir,'psclonlat.in');
if isfile(geoInput)
    lines=strtrim(splitlines(string(fileread(geoInput))));lines(lines=="")=[];
    if numel(lines)>=3
        for k=2:3
            file=strip(char(lines(k)),'both','"');[~,name,ext]=fileparts(file);
            copy(file,fullfile('export_geo',[name ext]));
        end
    end
end
if isfield(report,'processingAudit')&&isfield(report.processingAudit,'preprocessingFolder')
    processing=report.processingAudit.preprocessingFolder;
    for name={'CSK_subsetMgraph2run.xml','CSK_subsetgraph2run.xml','CSK_coreg_ifg2run.xml','CSK_export2run.xml'}
        copy(fullfile(processing,'graphs',name{1}),fullfile('graphs',name{1}));
    end
    exportRoot=fileparts(fileparts(lonPathFromInput(workDir)));
    copy(fullfile(exportRoot,'phase_preprocessing_aoi.json'),'phase_preprocessing_aoi.json');
    geometries=dir(fullfile(processing,'graphs','*.geometry.json'));
    for j=1:numel(geometries)
        copy(fullfile(geometries(j).folder,geometries(j).name),fullfile('graphs',geometries(j).name));
    end
    if isfield(report.config,'master_date')
        product(fullfile(processing,'master',[char(string(report.config.master_date)) '_sub.dim']),'master');
    end
    % The last export graph identifies the exact pair last written to geo/.
    graph=fullfile(processing,'graphs','CSK_export2run.xml');
    if isfile(graph)
        document=xmlread(graph);nodes=document.getElementsByTagName('node');
        for k=0:nodes.getLength()-1
            node=nodes.item(k);id=char(node.getAttribute('id'));
            if ~ismember(id,{'Read','Read(2)'}),continue;end
            files=node.getElementsByTagName('file');if files.getLength()==0,continue;end
            file=char(files.item(0).getTextContent());
            if strcmp(id,'Read'),label='coreg';else,label='interferogram';end
            product(strtrim(file),label);
        end
    end
    for name={'CSK_coreg_ifg_proc_stdout.log','CSK_export_proc_stdout.log','CSK_subset_proc_stdout.log','CSK_master_proc_stdout.log'}
        copy(fullfile(processing,'logs',name{1}),fullfile('logs',name{1}));
    end
end
destination=fullfile(fileparts(reportPath),'PHASE_CSK_diagnostics.zip');
zip(destination,{'*'},staging);

    function copy(source,relative)
        if ~isfile(source),return;end
        target=fullfile(staging,relative);parent=fileparts(target);
        if ~isfolder(parent),mkdir(parent);end
        [ok,message]=copyfile(source,target);
        if ~ok,error('PHASE:DiagnosticCopyFailed','Cannot collect %s: %s',source,message);end
    end
    function product(file,label)
        if ~isfile(file),return;end
        [parent,name]=fileparts(file);copy(file,fullfile(label,[name '.dim']));
        data=fullfile(parent,[name '.data']);
        grids=dir(fullfile(data,'tie_point_grids','*'));
        for j=1:numel(grids)
            if grids(j).isdir,continue;end
            copy(fullfile(grids(j).folder,grids(j).name),fullfile(label,[name '.data'],'tie_point_grids',grids(j).name));
        end
        if strcmp(label,'interferogram')
            for band={'orthorectifiedLat','orthorectifiedLon','elevation'}
                for ext={'.img','.hdr'}
                    copy(fullfile(data,[band{1} ext{1}]),fullfile(label,[name '.data'],[band{1} ext{1}]));
                end
            end
        end
    end
end

function file=lonPathFromInput(workDir)
lines=strtrim(splitlines(string(fileread(fullfile(workDir,'psclonlat.in')))));lines(lines=="")=[];
file=strip(char(lines(2)),'both','"');
end
