function cskRegressionSelfTest()
%CSKREGRESSIONSELFTEST Actual .m runner with controlled GPT and slave failure.
% This checks orchestration and logs; it is not a real-data SNAP validation.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'PHASE_Preprocessing'));
python=phase_preprocessing_beta.resolvePython('python3');
oldPath=path;folder=tempname;mkdir(folder);cleanup=onCleanup(@()cleanupTest(oldPath,folder)); %#ok<NASGU>
installed=fullfile(folder,'installed engine');bin=fullfile(installed,'PHASE_Preprocessing','snap2stamps','bin');mkdir(bin);
source=fullfile(root,'PHASE_Preprocessing','snap2stamps','bin');
for name={'CSK_master_selector.py','CSK_subset_master.py','phase_subprocess.py'}
    copyfile(fullfile(source,name{1}),fullfile(bin,name{1}));
end
selector=fullfile(bin,'CSK_master_selector.py');text=fileread(selector);
fid=fopen(selector,'w');fprintf(fid,'import io, json, urllib.request\nurllib.request.urlopen=lambda *a,**k:io.BytesIO(json.dumps({"daily":{"precipitation_sum":[0.0]}}).encode())\n%s',text);fclose(fid);
graphs=fullfile(installed,'PHASE_Preprocessing','snap2stamps','graphs');mkdir(graphs);
copyfile(fullfile(root,'PHASE_Preprocessing','snap2stamps','graphs','CSK_AOI_subset.xml'),graphs);
gpt=fullfile(folder,'fake_gpt');fid=fopen(gpt,'w');fprintf(fid,'#!%s\n%s',python,fileread(fullfile(root,'tests','fixtures','csk_gpt_stub.py')));fclose(fid);
fileattrib(gpt,'+x','a');
for name={'CSK_slaves_prep.py','CSK_subset_slaves.py','CSK_coreg_ifg.py','CSK_stamps_export.py','CSK_average_intensity.py','CSK_terrain_correction.py'}
    fid=fopen(fullfile(bin,name{1}),'w');fprintf(fid,'import sys\nprint("TEST SLAVE ENTRY REACHED",flush=True)\nsys.exit(77)\n');fclose(fid);
end
for automatic=[false true]
    project=fullfile(folder,sprintf('CSK run %d & data',automatic));mkdir(project);
    cfg=phase_preprocessing_beta.defaultConfig();cfg.constellation='CSK';cfg.python=python;
    cfg.auto_master=automatic;cfg.master_date='20210815';cfg.process_master=true;
    cfg.lon_min=9;cfg.lon_max=10;cfg.lat_min=45;cfg.lat_max=46;
    cfg.cpu=1;cfg.cache='128M';cfg.gptbin_path=gpt;
    phase_preprocessing_beta.saveConfig(project,cfg);
    slaves=fullfile(phase_preprocessing_beta.dataFolder(project),'slaves');mkdir(slaves);
    for date={'20210801','20210815','20210901'}
        fid=fopen(fullfile(slaves,sprintf('CSKS1_SCS_B_%s123456_%s123500.h5',date{1},date{1})),'w');fwrite(fid,'HDF fixture');fclose(fid);
    end
    engine=phase_preprocessing_beta.ProcessingEngine(project,installed);
    phase_preprocessing_beta.applyConfig(engine,cfg);caught=false;
    try,engine.run();catch e,caught=strcmp(e.identifier,'PHASE:ExternalCommandFailed');end
    assert(caught,'CSK runner failed before reaching the controlled slave operation.');
    assert(isfile(fullfile(phase_preprocessing_beta.dataFolder(project),'master','20210815_sub.dim')));
    log=fileread(engine.RunLogPath);
    assert(contains(log,'Master processing finished')&&contains(log,'TEST SLAVE ENTRY REACHED'));
    assert(contains(log,'PHASE:ExternalCommandFailed')&&~contains(log,'Unrecognized function or variable'));
end
% Even a failure loading the configuration must leave a persistent run log.
project=fullfile(folder,'missing config');mkdir(project);
engine=phase_preprocessing_beta.ProcessingEngine(project,installed);engine.constellation="CSK";
try,engine.run();catch,end
assert(isfile(engine.RunLogPath)&&contains(fileread(engine.RunLogPath),'input_preprocessing.mat'));
disp('CSK .m runner passed automatic/manual handoff and early-failure logging regressions.');
end

function cleanupTest(oldPath,folder)
path(oldPath);rmdir(folder,'s');
end
