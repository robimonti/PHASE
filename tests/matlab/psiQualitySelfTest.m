function psiQualitySelfTest()
%PSIQUALITYSELFTEST Freshness, portable paths, optional cleanup and editable exports.
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
export=fullfile(folder,'INSAR_20250101');mkdir(fullfile(export,'quality'));mkdir(fullfile(folder,'geo'));
source=fullfile(folder,'geo','test.lon');fid=fopen(source,'w');fwrite(fid,single(ones(4,1)),'single');fclose(fid);
info=dir(source);modified=double(java.io.File(source).lastModified())/1000;
sources=struct('path',source,'relativePath','geo\test.lon','bytes',info.bytes,'modified',modified,'requiredForFreshness',true);
sources(2)=struct('path',fullfile(folder,'removed.img'),'relativePath','removed.img','bytes',0,'modified',0,'requiredForFreshness',false);
report=struct('schema',1,'algorithm','psi-screening-v1','sources',sources, ...
    'grid',struct('bounds',[9 9.001 45 45.001],'rows',2,'columns',2), ...
    'acquisitionCount',6,'interferogramCount',5, ...
    'reasonCodes',zeros(2),'reasonBits',struct('bit',32,'message','Indicator unavailable'), ...
    'geometry',struct('status','unavailable','message','Insufficient terrain metadata.'));
for name={'suitability','geometry','stability','coverage','coherence'}
    values=[0 .3;.6 1];
    if strcmp(name{1},'suitability'),values=[0 1;2 3];end
    if strcmp(name{1},'geometry'),values(:)=nan;end
    report.layers.(name{1})=struct('title',name{1},'values',values);
end
path=fullfile(export,'quality','quality.json');fid=fopen(path,'w');fwrite(fid,jsonencode(report));fclose(fid);
data=phase_project.readPsiQuality(export);assert(data.available,data.message);
phase_preprocessing_beta.exportQualityFigure(data);
assert(isfile(fullfile(data.folder,'figures','PSI_quality.png')));
assert(isfile(fullfile(data.folder,'figures','PSI_quality.fig')));
invalid=report;invalid.layers.coverage.values(1)=2;
fid=fopen(path,'w');fwrite(fid,jsonencode(invalid));fclose(fid);
rejected=phase_project.readPsiQuality(export);assert(~rejected.available);
fid=fopen(path,'w');fwrite(fid,jsonencode(report));fclose(fid);
fid=fopen(source,'a');fwrite(fid,uint8(0));fclose(fid);
stale=phase_project.readPsiQuality(export);assert(~stale.available);
fprintf('PSI quality checks passed: portable paths, freshness, optional cleanup, PNG and FIG.\n');
end
