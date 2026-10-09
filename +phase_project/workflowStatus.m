function status = workflowStatus(projectRoot)
%WORKFLOWSTATUS Infer completed stages from published project products.
% A stage is green only when its expected output files exist. Opening a
% section or merely creating its working directory never marks it complete.

status = struct('preprocessing',false,'stamps',false,'model',false, ...
    'recommended','preprocessing');
if nargin < 1 || isempty(projectRoot), return; end
p = phase_project.paths(projectRoot);

datasets = [dir(fullfile(p.stamps,'ASC_*')); ...
    dir(fullfile(p.stamps,'DSC_*')); ...
    dir(fullfile(p.stamps,'DES_*'))];
datasets = datasets([datasets.isdir]);
exports = dir(fullfile(p.preprocessing,'INSAR_*'));
exports = exports([exports.isdir]);
for k = 1:numel(exports)
    folder = fullfile(exports(k).folder,exports(k).name);
    if ~isempty(datasets) && hasFiles(fullfile(folder,'diff0')) && ...
            hasFiles(fullfile(folder,'geo'))
        status.preprocessing = true;
        break
    end
end

published = dir(p.exports);
published = published([published.isdir] & ...
    ~ismember({published.name},{'.','..'}));
for k = 1:numel(published)
    folder = fullfile(published(k).folder,published(k).name);
    if ~isempty(dir(fullfile(folder,'*.csv'))) && ...
            ~isempty(dir(fullfile(folder,'*.xlsx')))
        status.stamps = true;
        break
    end
end

outputs = dir(fullfile(p.model,'output_*'));
outputs = outputs([outputs.isdir]);
for k = 1:numel(outputs)
    folder = fullfile(outputs(k).folder,outputs(k).name);
    if ~isempty(dir(fullfile(folder,'*.xlsx'))) && ...
            isfile(fullfile(folder,'files','mat','PHASEresults.mat'))
        status.model = true;
        break
    end
end

if ~status.preprocessing
    status.recommended = 'preprocessing';
elseif ~status.stamps
    status.recommended = 'stamps';
elseif ~status.model
    status.recommended = 'model';
else
    status.recommended = '';
end
end

function available = hasFiles(folder)
available = false;
if ~isfolder(folder), return; end
items = dir(folder);
available = any(~[items.isdir]);
end
