function stampsSafetySelfTest(stampsRoot)
%STAMPSSAFETYSELFTEST Exercise patched upstream files on temporary fixtures.
oldPath=path;oldDir=pwd;folder=tempname;mkdir(folder);
cleanup=onCleanup(@()finish(oldPath,oldDir,folder)); %#ok<NASGU>
runtime=fullfile(folder,'runtime');matlabDir=fullfile(runtime,'matlab');mkdir(matlabDir);
for name={'ps_select.m','ps_weed.m','ps_merge_patches.m','stamps.m'}
    copyfile(fullfile(stampsRoot,'matlab',name{1}),matlabDir);
end
r=phase_stamps_beta.repairStaMPSCompatibility(runtime);assert(r.changed);
first=fileread(fullfile(matlabDir,'stamps.m'));
r=phase_stamps_beta.repairStaMPSCompatibility(runtime);assert(~r.changed);
assert(strcmp(first,fileread(fullfile(matlabDir,'stamps.m'))));
assert(contains(first,"load('no_ps_info.mat','stamps_step_no_ps');"));
mock=fullfile(folder,'mocks');mkdir(mock);
writeText(fullfile(mock,'logit.m'),'function logit(varargin)\nend\n');
writeText(fullfile(mock,'getparm.m'),strjoin({ ...
    'function value=getparm(name,varargin)', ...
    'switch name', ...
    'case {''small_baseline_flag'',''weed_zero_elevation'',''weed_neighbours''},value=''n'';', ...
    'case {''merge_resample_size''},value=0;', ...
    'case {''merge_standard_dev'',''weed_max_noise''},value=Inf;', ...
    'case ''drop_ifg_index'',value=[];', ...
    'case ''weed_time_win'',value=730;', ...
    'case ''weed_standard_dev'',value=1;', ...
    'otherwise,error(''Unexpected mocked parameter: %s'',name);', ...
    'end','end'},newline));
addpath(matlabDir,'-begin');addpath(mock,'-begin');
work=fullfile(folder,'work');mkdir(work);cd(work);
% Empty patch.list previously produced a misleading merge/origin exception.
writeText('patch.list','');
mustFail(@()ps_merge_patches(),'PHASE_StaMPS:noPersistentScatterers');
% Actual upstream ps_weed reaches the guard with two selected points.
psver=1;save('psver.mat','psver');
ps=struct('n_ifg',2,'n_ps',2,'day',[1;2],'bperp',[0;10], ...
    'master_day',1,'xy',[1 0 0;2 10 10],'ij',[1 0 0;2 1 1], ...
    'lonlat',[9 45;9.001 45.001],'sort_ix',[1;2]);
save('ps1.mat','-struct','ps');
ph=ones(2,2);save('ph1.mat','ph');ph_patch=ones(2,2);save('pm1.mat','ph_patch');
ix=[1;2];keep_ix=true(2,1);K_ps2=zeros(2,1);C_ps2=K_ps2;coh_ps2=ones(2,1);
save('select1.mat','ix','keep_ix','K_ps2','C_ps2','coh_ps2');
mustFail(@()ps_weed(0,1,0),'PHASE_StaMPS:insufficientPersistentScatterers');
% Verify the caller no longer overwrites the flag saved by a called step.
clear stamps ps_weed
mini=fullfile(folder,'mini','matlab');mkdir(mini);
writeText(fullfile(mini,'stamps.m'),strjoin({ ...
    'function stamps()', ...
    'stamps_step_no_ps=zeros(5,1);', ...
    'ps_weed;', ...
    'save(''no_ps_info.mat'',''stamps_step_no_ps'');', ...
    'end'},newline));
writeText(fullfile(mini,'ps_weed.m'),strjoin({ ...
    'function ps_weed()', ...
    'stamps_step_no_ps=zeros(5,1);stamps_step_no_ps(4)=1;', ...
    'save(''no_ps_info.mat'',''stamps_step_no_ps'');', ...
    'end'},newline));
phase_stamps_beta.repairStaMPSCompatibility(fileparts(mini));
addpath(mini,'-begin');stamps();flags=load('no_ps_info.mat');
assert(flags.stamps_step_no_ps(4)==1,'The caller overwrote the no-PS flag.');
fprintf('StaMPS safety regressions passed: empty merge, two-point weeding, saved flags and idempotency.\n');
end
function writeText(file,text)
fid=fopen(file,'w');fprintf(fid,'%s',strrep(text,'\n',newline));fclose(fid);
end
function mustFail(action,id)
try,action();catch ME,assert(strcmp(ME.identifier,id),'Unexpected error [%s]: %s',ME.identifier,ME.message);return;end
error('Expected failure: %s',id);
end
function finish(oldPath,oldDir,folder)
cd(oldDir);path(oldPath);clear stamps ps_weed ps_merge_patches getparm logit;rmdir(folder,'s');
end
