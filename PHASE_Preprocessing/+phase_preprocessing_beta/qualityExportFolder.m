function [folder,master] = qualityExportFolder(projectDir,cfg)
%QUALITYEXPORTFOLDER Use the master actually passed to SNAP, including auto selection.
project=phase_preprocessing_beta.dataFolder(projectDir);master=char(string(cfg.master_date));
conf=fullfile(project,'snap2stamps','bin','project.conf');
if isfile(conf)
    token=regexp(fileread(conf),'(?m)^\s*MASTER\s*=\s*(\d{8})\s*$','tokens','once');
    if ~isempty(token),master=token{1};end
end
folder=fullfile(project,['INSAR_' master]);
end
