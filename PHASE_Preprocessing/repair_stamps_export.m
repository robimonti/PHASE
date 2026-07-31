function outputs = repair_stamps_export(workDir)
%REPAIR_STAMPS_EXPORT Rebuild StaMPS time-series exports without truncation.
% Uses completed ps_plot MAT products; it does not rerun StaMPS processing.
if nargin < 1 || isempty(workDir)
    workDir = uigetdir(pwd, 'Select the StaMPS processing folder');
    if isequal(workDir, 0), error('PHASE:repairExportCancelled', 'No folder selected.'); end
end
workDir = char(workDir);
cfg = load(fullfile(workDir, 'input_StaMPS.mat'), 'export_name', 'year_0', 'month_0', 'day_0');
exportDir = fullfile(workDir, 'EXPORT');
if ~isfolder(exportDir), mkdir(exportDir); end

types = {'v-dao', 'v-do'};
tsPath = ''; velPath = '';
for k = 1:numel(types)
    candidateTs = fullfile(workDir, ['ps_plot_ts_' types{k} '.mat']);
    candidateVel = fullfile(workDir, ['ps_plot_' types{k} '.mat']);
    if exist(candidateTs, 'file') == 2 && exist(candidateVel, 'file') == 2
        tsPath = candidateTs; velPath = candidateVel; valueType = types{k}; break
    end
end
if isempty(tsPath)
    error('PHASE:repairExportProductsMissing', ...
        'Matching ps_plot_ts_* and ps_plot_v-* MAT files were not found in %s.', workDir);
end
ts = load(tsPath); vel = load(velPath);
assertFields(ts, {'day', 'lonlat', 'master_day', 'ph_mm'}, tsPath);
assertFields(vel, {'ph_disp'}, velPath);

timeZero = datenum(cfg.year_0, cfg.month_0, cfg.day_0);
timeDays = double(ts.day(:)') - timeZero;
timeMaster = double(ts.master_day) - timeZero;
dispTs = ts.ph_mm;
if size(dispTs, 2) ~= numel(timeDays)
    error('PHASE:repairExportShape', 'ph_mm and day have inconsistent column counts.');
end
beforeMaster = timeDays < timeMaster;
timeDays = [timeDays(beforeMaster), timeMaster, timeDays(~beforeMaster)];
dispTs = [dispTs(:, beforeMaster), zeros(size(dispTs, 1), 1), dispTs(:, ~beforeMaster)];
if size(dispTs, 1) ~= size(ts.lonlat, 1) || size(vel.ph_disp, 1) ~= size(ts.lonlat, 1)
    error('PHASE:repairExportShape', 'The MAT products have inconsistent PS counts.');
end

base = fullfile(exportDir, char(string(cfg.export_name)));
writeSeries([ts.lonlat(:, 1), ts.lonlat(:, 2), vel.ph_disp, dispTs], ...
    timeDays, [base '.xlsx'], [base '.csv']);
ind = true(size(ts.lonlat, 1), 1);
save(fullfile(workDir, 'PS_index.mat'), 'ind');
outputs = {[base '.xlsx'], [base '.csv']};
fprintf('Repaired %s export: %d PS, %d time columns.\n', ...
    valueType, size(dispTs, 1), size(dispTs, 2));
end

function writeSeries(data, timeDays, xlsxPath, csvPath)
if exist(xlsxPath, 'file') == 2, delete(xlsxPath); end
point = (1:size(data, 1))';
infos = table(point, data(:, 1), data(:, 2), data(:, 3), ...
    'VariableNames', {'point', 'lon', 'lat', 'vel'});
writetable(infos, xlsxPath, 'Range', 'A2');
% Start-cell-only range writes every supplied column; there is no FZ limit.
series = table([timeDays; data(:, 4:end)]);
writetable(series, xlsxPath, 'Range', 'E1');
writetable(readtable(xlsxPath), csvPath, 'WriteMode', 'overwrite');
end

function assertFields(s, fields, pathName)
for k = 1:numel(fields)
    if ~isfield(s, fields{k})
        error('PHASE:repairExportMissingField', 'Missing %s in %s.', fields{k}, pathName);
    end
end
end
