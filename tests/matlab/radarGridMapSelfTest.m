function radarGridMapSelfTest(outputDir)
%RADARGRIDMAPSELFTEST True cell boundaries, stages, bounded viewport and exports.
root=tempname;mkdir(root);cleanup=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
patch=fullfile(root,'PATCH_1');mkdir(patch);
width=4;len=3;write(fullfile(root,'width.txt'),'4');write(fullfile(root,'len.txt'),'3');
[rg,az]=meshgrid(0:width-1,0:len-1);
lon=9+rg*.001;lat=45+az*.001;
for item={{'lon',lon},{'lat',lat}}
    pair=item{1};file=fullfile(root,['geo.' pair{1}]);fid=fopen(file,'w','ieee-be');
    fwrite(fid,pair{2}','single');fclose(fid);
end
write(fullfile(root,'psclonlat.in'),sprintf('4\n%s\n%s\n',fullfile(root,'geo.lon'),fullfile(root,'geo.lat')));
ij=[1 0 0;2 1 1;3 2 3];lonlat=[9 45;9.001 45.001;9.003 45.002];n_ps=3;
save(fullfile(patch,'ps1.mat'),'ij','lonlat','n_ps');
D_A=[.1;.2;.3];save(fullfile(patch,'da1.mat'),'D_A');
ix=[1;2;3];keep_ix=[false;true;true];save(fullfile(patch,'select1.mat'),'ix','keep_ix');
ix_weed=[false;true];save(fullfile(patch,'weed1.mat'),'ix_weed');
cfg=phase_stamps_beta.defaultConfig();cfg.export_name='fixture';
legacy=cfg;legacy.export_grid_diagnostics=true;
save(fullfile(root,'input_StaMPS.mat'),'-struct','legacy');
loaded=phase_stamps_beta.loadConfig(root);
assert(~isfield(loaded,'export_grid_diagnostics'));
assert(strcmp(loaded.export_name,cfg.export_name));
source=phase_stamps_beta.exportGridDiagnostics(root,cfg,[],struct('interactiveOnly',true));
assert(~isfolder(fullfile(root,'EXPORT')),'Viewing must not write diagnostic exports.');
assert(isequal(source.candidates.stage,[1;2;3]));
frame=phase_stamps_beta.radarGridMap(source,struct('bounds',[8.99 9.01 44.99 45.01]));
assert(frame.visibleCount==3&&isequal(frame.stageCounts,[1 1 1]));
assert(frame.azimuthStep==1&&frame.rangeStep==1);
% Grid lines run at half-pixel boundaries around zero-based pixel centres.
assert(max(abs(frame.lines{1}(1,:)-[8.9995 44.9995]))<2e-6);
% Missing geolocation breaks curves rather than drawing a line to (0,0).
withHole=source;withHole.geo.lon(3,2)=NaN;withHole.geo.lat(3,2)=NaN;
withHole.geo.lonF.Values=withHole.geo.lon;withHole.geo.latF.Values=withHole.geo.lat;
hole=phase_stamps_beta.radarGridMap(withHole,struct('bounds',[8.99 9.01 44.99 45.01]));
assert(any(cellfun(@(line)any(isnan(line),'all'),hole.lines)));
one=phase_stamps_beta.radarGridMap(source,struct('stage',3));assert(size(one.points,1)==1&&one.points(1,4)==3);
empty=phase_stamps_beta.radarGridMap(source,struct('bounds',[10 11 46 47]));
assert(empty.visibleCount==0&&isempty(empty.lines));
% The browser payload stays bounded on large candidate sets without changing counts.
names=fieldnames(source.candidates);
for k=1:numel(names),source.candidates.(names{k})=repmat(source.candidates.(names{k}),15000,1);end
many=phase_stamps_beta.radarGridMap(source,struct('bounds',[8.99 9.01 44.99 45.01]));
assert(many.visibleCount==45000&&many.displayedCount==40000&&size(many.points,1)==40000);
try
    phase_stamps_beta.radarGridMap(source,struct('bounds',[11 10 46 47]));
    error('Test:ExpectedFailure','Invalid bounds were accepted.');
catch ME
    assert(strcmp(ME.identifier,'PHASE:GridMapBounds'));
end
if nargin>0
    write(fullfile(outputDir,'phase-radar-grid-ui-fixture.json'),jsonencode(frame));
    exported=phase_stamps_beta.exportGridDiagnostics(root,cfg,[],struct('bounds',frame.bounds));
    assert(isfile(exported.daFigure)&&isfile(exported.stageFigure));
    assert(isfile(strrep(exported.daFigure,'.png','.fig')));
    assert(isfile(strrep(exported.stageFigure,'.png','.fig')));
    allCandidates=readtable(exported.csv);assert(height(allCandidates)==3);
end
fprintf('Radar-grid source, cell boundaries, stage counts, empty views and bounded large datasets passed.\n');
end

function write(path,text)
fid=fopen(path,'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',text);
end
