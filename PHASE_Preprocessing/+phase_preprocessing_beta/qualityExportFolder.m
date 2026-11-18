function [folder,master,info] = qualityExportFolder(projectDir,cfg)
%QUALITYEXPORTFOLDER Resolve the actual SNAP master without choosing the newest dataset.
% MASTER normally contains a processed product path, not a bare date. Parse
% only its basename so dates in parent folders cannot change the reference.
project=phase_preprocessing_beta.dataFolder(projectDir);
master=strtrim(char(string(cfg.master_date)));folder=fullfile(project,['INSAR_' master]);
info=struct('source','saved-config','message','Using the saved master date.');
conf=fullfile(project,'snap2stamps','bin','project.conf');
if isfile(conf)
    token=regexp(fileread(conf),'(?m)^\s*MASTER\s*=\s*([^\r\n]+)','tokens','once');
    if ~isempty(token)
        actual=masterDate(token{1});
        if ~isempty(actual)
            master=actual;folder=fullfile(project,['INSAR_' master]);
            info=struct('source','project.conf','message','Using the master product actually passed to SNAP.');return;
        end
    end
end
% Older projects may retain the processed master without project.conf.
products=dir(fullfile(project,'master','*.dim'));dates={};
for k=1:numel(products)
    actual=masterDate(products(k).name);
    if ~isempty(actual),dates{end+1}=actual;end %#ok<AGROW>
end
dates=unique(dates);
if numel(dates)==1 && completeExport(project,dates{1})
    master=dates{1};folder=fullfile(project,['INSAR_' master]);
    info=struct('source','master-product','message','Using the retained processed master product.');return;
end
manual=~isfield(cfg,'auto_master')||~logical(cfg.auto_master);
if manual && completeExport(project,master),return;end
exports=dir(fullfile(project,'INSAR_*'));complete={};
for k=1:numel(exports)
    date=regexp(exports(k).name,'^INSAR_(\d{8})$','tokens','once');
    if exports(k).isdir && ~isempty(date) && completeExport(project,date{1})
        complete{end+1}=date{1}; %#ok<AGROW>
    end
end
if numel(complete)==1
    master=complete{1};folder=fullfile(project,['INSAR_' master]);
    info=struct('source','single-export','message','Using the only complete SNAP export in this project.');
elseif numel(complete)>1
    folder='';master='';
    info=struct('source','ambiguous','message', ...
        ['More than one complete SNAP export is available (' strjoin(sort(complete),', ') ...
        '). Choose the intended master manually in Preprocessing and save, or restore its project.conf before building the quality map.']);
end
end

function date=masterDate(value)
value=strtrim(char(string(value)));value=strrep(value,'\','/');
value=regexprep(value,'^["'']|["'']$','');parts=strsplit(value,'/');name=parts{end};
date='';
patterns={'^(\d{8})$','^(\d{8})(?:_sub)?\.dim$', ...
    '_(\d{8})T\d{6}(?:_|\.)','_(\d{8})\d{6}(?:_|\.)'};
for k=1:numel(patterns)
    token=regexp(name,patterns{k},'tokens','once');
    if ~isempty(token),date=token{1};return;end
end
end

function valid=completeExport(project,date)
valid=false;
if isempty(regexp(date,'^\d{8}$','once')),return;end
root=fullfile(project,['INSAR_' date]);
pars=dir(fullfile(root,'*slc',[date '.*slc.par']));
lon=dir(fullfile(root,'geo','*.lon'));lat=dir(fullfile(root,'geo','*.lat'));
valid=~isempty(pars)&&~isempty(lon)&&~isempty(lat);
end
