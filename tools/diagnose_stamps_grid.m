function report = diagnose_stamps_grid(workDir,outputPath)
%DIAGNOSE_STAMPS_GRID Collect small, read-only PSI/geolocation diagnostics.
% Select a PSI dataset containing input_StaMPS.mat and/or PATCH_*.
% Writes a JSON/ZIP diagnostic report only; no inputs, parameters or results change.
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
report=struct('diagnosticSchema',3,'generatedAt',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
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
try
    report.coverageAudit=coverageAudit(workDir,report);
catch ME
    report.coverageAuditError=ME.message;
end
[selection,selectionErrors]=selectionStages(workDir);
report.selectionStages=cellfun(@(v)struct('stage',v.name,'points',size(v.coordinates,1),'modified',v.modified),selection,'UniformOutput',false);
report.selectionStageErrors=selectionErrors;
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
[stages,~]=selectionStages(workDir);
if ~isempty(stages),mkdir(fullfile(staging,'selection_stages'));end
for j=1:numel(stages)
    item=stages{j};data=array2table([item.coordinates item.metrics], ...
        'VariableNames',[{'longitude','latitude'} item.metricNames]);
    writetable(data,fullfile(staging,'selection_stages',[item.name '.csv']));
end
for name={'psi_area.json','psi_geolocation.json','psclonlat.in','selpsc.in','calamp.in','calamp.out','patch.list'}
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
    input=fullfile(processing,'input_preprocessing.mat');
    if isfile(input)
        vars=whos('-file',input);
        if ismember('dem_file',{vars.name})
            cfg=load(input,'dem_file');info=dir(cfg.dem_file);
            if ~isempty(info)&&info.bytes<=100*1024^2
                copy(cfg.dem_file,'dem/source.tif');
            end
        end
    end
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
            coherence=dir(fullfile(data,'coh*'));
            for j=1:numel(coherence)
                [~,~,extension]=fileparts(coherence(j).name);
                if ~coherence(j).isdir&&ismember(extension,{'.img','.hdr'})&&coherence(j).bytes<=100*1024^2
                    copy(fullfile(data,coherence(j).name),fullfile(label,[name '.data'],coherence(j).name));
                end
            end
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

function [stages,errors]=selectionStages(workDir)
% Load coordinates/quality vectors only, never a phase or time-series matrix.
stages={};errors={};patches=dir(fullfile(workDir,'PATCH_*'));patches=patches([patches.isdir]);
folders=[{workDir},arrayfun(@(p)fullfile(p.folder,p.name),patches,'UniformOutput',false)];
labels=[{'ROOT'},{patches.name}];
for k=1:numel(folders)
    folder=folders{k};label=labels{k};
    try
        initial=loadSmall(fullfile(folder,'ps1.mat'),{'lonlat'});
        if isfield(initial,'lonlat')
            xy=double(initial.lonlat);n=size(xy,1);
            da=loadSmall(fullfile(folder,'da1.mat'),{'D_A'});
            pm=loadSmall(fullfile(folder,'pm1.mat'),{'coh_ps'});
            add('candidates',xy,{'amplitude_dispersion','initial_coherence'},[vector(da,'D_A',n) vector(pm,'coh_ps',n)],'ps1.mat');
            sl=loadSmall(fullfile(folder,'select1.mat'),{'ix','ix2','keep_ix','coh_ps2','coh_thresh'});
            if isfield(sl,'ix')
                ix=indices(sl.ix,n);count=numel(ix);
                quality=[vector(sl,'coh_ps2',count) vector(sl,'coh_thresh',count)];
                add('selected',xy(ix,:),{'reestimated_coherence','coherence_threshold'},quality,'select1.mat');
                if isfield(sl,'keep_ix')
                    keep=sl.keep_ix(:);
                    if numel(keep)~=count||any(~ismember(keep,[0 1])),error('Invalid keep_ix size/values.');end
                    keep=logical(keep);
                    add('reselected',xy(ix(keep),:),{'reestimated_coherence','coherence_threshold'},quality(keep,:),'select1.mat');
                end
            elseif isfield(sl,'ix2')
                ix=indices(sl.ix2,n);
                add('reselected',xy(ix,:),{'reestimated_coherence','coherence_threshold'},[vector(sl,'coh_ps2',numel(ix)) vector(sl,'coh_thresh',numel(ix))],'select1.mat');
            end
        end
        final=loadSmall(fullfile(folder,'ps2.mat'),{'lonlat'});
        if isfield(final,'lonlat'),add('retained',double(final.lonlat),{},zeros(size(final.lonlat,1),0),'ps2.mat');end
        if k==1
            files=dir(fullfile(folder,'ps_plot_ts_v*.mat'));
            for j=1:numel(files)
                exported=loadSmall(fullfile(files(j).folder,files(j).name),{'lonlat'});
                if isfield(exported,'lonlat')
                    [~,stem]=fileparts(files(j).name);
                    add(['export_' stem],double(exported.lonlat),{},zeros(size(exported.lonlat,1),0),files(j).name);
                end
            end
            files=dir(fullfile(folder,'EXPORT','*.csv'));
            for j=1:numel(files)
                file=fullfile(files(j).folder,files(j).name);opts=detectImportOptions(file,'VariableNamingRule','preserve');
                names=lower(string(opts.VariableNames));
                ilon=find(ismember(names,["lon","longitude"]),1);ilat=find(ismember(names,["lat","latitude"]),1);
                if isempty(ilon)||isempty(ilat)
                    if numel(names)<3,continue;end
                    ilon=2;ilat=3; % PHASE exports: point, longitude, latitude, ...
                end
                opts.SelectedVariableNames=opts.VariableNames([ilon ilat]);tableData=readtable(file,opts);
                xy=nan(height(tableData),2);
                for axis=1:2
                    v=tableData{:,axis};if isnumeric(v),xy(:,axis)=double(v);else,xy(:,axis)=str2double(string(v));end
                end
                valid=all(isfinite(xy),2)&abs(xy(:,1))<=180&abs(xy(:,2))<=90;
                xy=xy(valid,:);
                [~,stem]=fileparts(file);add(['csv_' stem],xy,{},zeros(size(xy,1),0),fullfile('EXPORT',files(j).name));
            end
        end
    catch ME
        errors{end+1}=struct('folder',folder,'error',ME.message); %#ok<AGROW>
    end
end
    function add(name,coordinates,metricNames,metrics,source)
        if size(coordinates,2)~=2,error('Invalid lonlat dimensions in %s.',source);end
        stage=struct('name',regexprep([label '_' name],'[^a-zA-Z0-9_]','_'),'coordinates',coordinates, ...
            'metricNames',{metricNames},'metrics',metrics,'modified',dirDate(fullfile(folder,source)));
        stages{end+1}=stage;
    end
end

function values=vector(data,name,n)
values=nan(n,1);
if ~isfield(data,name),return;end
v=double(data.(name)(:));
if isscalar(v),values(:)=v;elseif numel(v)==n,values=v;end
end

function ix=indices(values,n)
if islogical(values)
    if numel(values)~=n,error('Candidate selection mask length does not match ps1.');end
    ix=find(values);
else
    ix=double(values(:));
    if any(~isfinite(ix)|ix<1|ix>n|ix~=floor(ix)),error('Invalid candidate selection indices.');end
end
end

function result=loadSmall(file,names)
result=struct();if ~isfile(file),return;end
vars=whos('-file',file);fields=intersect(names,{vars.name},'stable');
if ~isempty(fields),result=load(file,fields{:});end
end

function out=coverageAudit(workDir,report)
% Compare the actual input mask to the current expected rule and inspect
% common nonzero/finite SLC support. No amplitudes/phases are exported.
out=struct('gridOrder','Rows: south to north; columns: west to east. Counts are radar pixels, not ground-area fractions.');
input=fullfile(workDir,'psclonlat.in');if ~isfile(input),out.message='No psclonlat.in retained.';return;end
lines=strtrim(splitlines(string(fileread(input))));lines(lines=="")=[];
if numel(lines)<3,out.message='Incomplete psclonlat.in.';return;end
lonPath=strip(char(lines(2)),'both','"');latPath=strip(char(lines(3)),'both','"');
if ~isfile(lonPath)||~isfile(latPath),out.message='The original geolocation rasters are unavailable.';return;end
exportRoot=fileparts(fileparts(lonPath));area=struct();areaFile=fullfile(exportRoot,'phase_preprocessing_aoi.json');
if isfile(areaFile),area=jsondecode(fileread(areaFile));end
if isfield(area,'requestedBounds')
    bounds=double(area.requestedBounds(:).');
elseif isfield(report,'processingAudit')&&isfield(report.processingAudit,'preprocessingConfig')&& ...
        all(isfield(report.processingAudit.preprocessingConfig,{'lon_min','lon_max','lat_min','lat_max'}))
    c=report.processingAudit.preprocessingConfig;bounds=[c.lon_min c.lon_max c.lat_min c.lat_max];
else,out.message='No geographic bounds available for spatial diagnostics.';return;end
if any(~isfinite(bounds))||bounds(1)>=bounds(2)||bounds(3)>=bounds(4),error('Invalid diagnostic geographic bounds.');end
out.bounds=bounds;nBins=64;out.bins=nBins;
n=floor(dirBytes(lonPath)/4);if dirBytes(latPath)~=n*4,error('Geolocation sizes disagree.');end
binIndex=zeros(n,1,'uint16');sourceCounts=zeros(nBins^2,1);includedCounts=sourceCounts;
maskName='psi_geolocation.mask';
if isfield(report.config,'psi_area_mode')&&strcmp(report.config.psi_area_mode,'subarea'),maskName='psi_area.mask';end
maskPath=fullfile(workDir,maskName);hasMask=isfile(maskPath);
if hasMask&&dirBytes(maskPath)~=n,error('The current PSI mask has a different size from the export.');end
out.maskAudit=struct('maskPresent',hasMask,'file',maskName,'unexpectedlyExcludedPixels',0,'unexpectedlyIncludedPixels',0, ...
    'note','Compares the retained mask to the current saved configuration, not a replay of historical preparation.');
lo=fopen(lonPath,'r','ieee-be');la=fopen(latPath,'r','ieee-be');mask=-1;
if hasMask,mask=fopen(maskPath,'r');end
cleanup=onCleanup(@()closeMany([lo la mask])); %#ok<NASGU>
if lo<0||la<0||(hasMask&&mask<0),error('Cannot open geolocation/mask diagnostics.');end
polygon=[];
if isfield(report.config,'psi_area_mode')&&strcmp(report.config.psi_area_mode,'subarea')&&isfield(report.config,'psi_polygon')
    polygon=double(report.config.psi_polygon);
    if size(polygon,2)~=2||size(polygon,1)<3,error('Invalid saved PSI polygon.');end
    polygon(:,1)=rad2deg(unwrap(deg2rad(polygon(:,1))));centre=mean(polygon(:,1));
end
for start=1:500000:n
    count=min(500000,n-start+1);x=fread(lo,count,'single=>double');y=fread(la,count,'single=>double');
    valid=isfinite(x)&isfinite(y)&abs(x)<=180&abs(y)<90&~(x==0&y==0);
    inBox=valid&x>=bounds(1)&x<=bounds(2)&y>=bounds(3)&y<=bounds(4);
    index=zeros(count,1,'uint16');
    col=min(nBins,floor((x(inBox)-bounds(1))/(bounds(2)-bounds(1))*nBins)+1);
    row=min(nBins,floor((y(inBox)-bounds(3))/(bounds(4)-bounds(3))*nBins)+1);
    index(inBox)=uint16(sub2ind([nBins nBins],row,col));binIndex(start:start+count-1)=index;
    sourceCounts=sourceCounts+accumarray(double(index(index>0)),1,[nBins^2 1]);
    expected=valid;
    if isfield(area,'requestedBounds'),expected=expected&inBox;end
    if ~isempty(polygon)
        shifted=centre+mod(x-centre+180,360)-180;
        expected=expected&inpolygon(shifted,y,polygon(:,1),polygon(:,2));
    end
    if hasMask
        actual=fread(mask,count,'uint8')==0;
        if numel(actual)~=count,error('Incomplete mask.');end
    else,actual=true(count,1);end
    out.maskAudit.unexpectedlyExcludedPixels=out.maskAudit.unexpectedlyExcludedPixels+nnz(expected&~actual);
    out.maskAudit.unexpectedlyIncludedPixels=out.maskAudit.unexpectedlyIncludedPixels+nnz(~expected&actual);
    selected=index>0&actual;
    includedCounts=includedCounts+accumarray(double(index(selected)),1,[nBins^2 1]);
end
clear cleanup
out.validGeolocationCounts=reshape(sourceCounts,nBins,nBins);
out.currentMaskIncludedCounts=reshape(includedCounts,nBins,nBins);
[stages,~]=selectionStages(workDir);out.selectionCounts={};
for k=1:numel(stages)
    p=stages{k}.coordinates;valid=all(isfinite(p),2);
    counts=histcounts2(p(valid,2),p(valid,1),linspace(bounds(3),bounds(4),nBins+1),linspace(bounds(1),bounds(2),nBins+1));
    out.selectionCounts{end+1}=struct('stage',stages{k}.name,'counts',counts); %#ok<AGROW>
end
calamp=fullfile(workDir,'calamp.in');out.slcSupport={};
if ~isfile(calamp),out.slcSupportNote='calamp.in not retained; common stack support could not be audited.';return;end
paths=strtrim(splitlines(string(fileread(calamp))));paths(paths=="")=[];
common=true(n,1);complete=true;
for k=1:numel(paths)
    file=strip(char(paths(k)),'both','"');if ~isfile(file),file=fullfile(workDir,file);end
    entry=struct('file',char(paths(k)),'status','unavailable');
    if ~isfile(file)||dirBytes(file)~=8*n
        complete=false;entry.status='Missing input or unsupported size (expected SNAP complex float32).';out.slcSupport{end+1}=entry;continue; %#ok<AGROW>
    end
    fprintf('Inspecting retained SLC support %d/%d (read only)...\n',k,numel(paths));
    fid=fopen(file,'r','ieee-be');if fid<0,error('Cannot read %s.',file);end
    cleanup=onCleanup(@()fclose(fid));validCount=0;counts=zeros(nBins^2,1);
    for start=1:500000:n
        count=min(500000,n-start+1);v=fread(fid,[2 count],'single=>double');
        if size(v,2)~=count,error('Incomplete complex SLC.');end
        valid=all(isfinite(v),1).'&any(v~=0,1).';validCount=validCount+nnz(valid);
        indices=start:start+count-1;common(indices)=common(indices)&valid;
        bins=binIndex(indices);keep=bins>0&valid;
        counts=counts+accumarray(double(bins(keep)),1,[nBins^2 1]);
    end
    clear cleanup
    entry.status='read';entry.finiteNonzeroPixels=validCount;entry.counts=reshape(counts,nBins,nBins);out.slcSupport{end+1}=entry; %#ok<AGROW>
end
if complete&&~isempty(paths)
    counts=accumarray(double(binIndex(common&binIndex>0)),1,[nBins^2 1]);
    out.commonSlcSupportCounts=reshape(counts,nBins,nBins);
    out.commonSlcSupportNote='Finite, nonzero complex samples in every listed acquisition. This does not establish phase coherence or PS suitability.';
else,out.slcSupportNote='Incomplete source stack: common support is unavailable; omitted files are not assumed valid.';end
end

function closeMany(handles)
for handle=handles
    if handle>=0,try,fclose(handle);catch,end,end
end
end
