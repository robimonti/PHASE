function psiQualityExportSelfTest()
%PSIQUALITYEXPORTSELFTEST Real MASTER path formats and conservative fallbacks.
root=tempname;mkdir(root);cleanup=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
project=fullfile(root,'PHASE_Preprocessing');mkdir(project);
conf=fullfile(project,'snap2stamps','bin','project.conf');mkdir(fileparts(conf));
cfg=phase_preprocessing_beta.defaultConfig();cfg.master_date='20200722';cfg.auto_master=true;
makeExport(project,'20251011');
paths={'D:\Tesi\20240101\master\S1A_IW_SLC__1SDV_20251011T045307_20251011T045334_split_IW1_Orb.dim', ...
    '"D:\Project with spaces\master\20251011_sub.dim"', ...
    '/mnt/data/master/S1B_IW_SLC__1SDV_20251011T045307_20251011T045334_split_IW2_Orb.dim', ...
    '/data/CSKS4_SCS_B_HI_01_HH_RA_SF_20251011045307_20251011045313.h5', ...
    '20251011'};
for k=1:numel(paths)
    fid=fopen(conf,'w');fprintf(fid,'PROJECTFOLDER=%s\r\nMASTER=%s\r\nCPU=4\r\n',project,paths{k});fclose(fid);
    [folder,master,info]=phase_preprocessing_beta.qualityExportFolder(root,cfg);
    assert(strcmp(master,'20251011')&&strcmp(info.source,'project.conf'));
    assert(strcmp(folder,fullfile(project,'INSAR_20251011')));
end
% A second older export must never override an explicit SNAP reference.
makeExport(project,'20200722');
[~,master]=phase_preprocessing_beta.qualityExportFolder(root,cfg);assert(strcmp(master,'20251011'));
delete(conf);
[folder,master,info]=phase_preprocessing_beta.qualityExportFolder(root,cfg);
assert(isempty(folder)&&isempty(master)&&strcmp(info.source,'ambiguous'));
cfg.auto_master=false;
[~,master,info]=phase_preprocessing_beta.qualityExportFolder(root,cfg);
assert(strcmp(master,'20200722')&&strcmp(info.source,'saved-config'));
cfg.auto_master=true;
rmdir(fullfile(project,'INSAR_20200722'),'s');
[~,master,info]=phase_preprocessing_beta.qualityExportFolder(root,cfg);
assert(strcmp(master,'20251011')&&strcmp(info.source,'single-export'));
mkdir(fullfile(project,'master'));fid=fopen(fullfile(project,'master','20251011_sub.dim'),'w');fclose(fid);
[~,master,info]=phase_preprocessing_beta.qualityExportFolder(root,cfg);
assert(strcmp(master,'20251011')&&strcmp(info.source,'master-product'));
% An explicit but incomplete new reference is not replaced by an old export.
fid=fopen(conf,'w');fprintf(fid,'MASTER=D:\\master\\20260912_sub.dim\r\n');fclose(fid);
[folder,master]=phase_preprocessing_beta.qualityExportFolder(root,cfg);
assert(strcmp(master,'20260912')&&strcmp(folder,fullfile(project,'INSAR_20260912')));
fprintf('Quality export resolution passed: S1/CSK paths, quotes, stale dates, manual selection and ambiguous exports.\n');
end

function makeExport(project,date)
folder=fullfile(project,['INSAR_' date]);mkdir(fullfile(folder,'rslc'));mkdir(fullfile(folder,'geo'));
fid=fopen(fullfile(folder,'rslc',[date '.rslc.par']),'w');fprintf(fid,'range_samples: 1\nazimuth_lines: 1\n');fclose(fid);
for suffix={'lon','lat'}
    fid=fopen(fullfile(folder,'geo',[date '.' suffix{1}]),'w','ieee-be');fwrite(fid,single(1),'single');fclose(fid);
end
end
