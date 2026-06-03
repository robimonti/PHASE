% launch_phase_with_demo_workdir.m
%
% Boots PHASE_StaMPS pre-configured with a synthetic project folder
% containing a fake ps_plot_ts_v-do.mat, so the new "TS Points" tab
% can be exercised end-to-end without a real PSI workdir.
%
% After the app opens, the active tab is switched to "TS Points" and
% the user only has to click the green "Open TS Points picker..."
% button.

phase_root  = fileparts(fileparts(mfilename('fullpath')));
stamps_root = fullfile(phase_root, 'StaMPS');

addpath(fullfile(stamps_root, 'matlab'));
addpath(fullfile(stamps_root, 'matlab_compat'));

% Synthetic workdir with a realistic ps_plot_ts_v-do.mat
tmp = tempname;
mkdir(tmp);
rng(42);

n_cluster = 200;
n_sparse  = 400;
lonlat = [16.60 + 0.005*randn(n_cluster,1), ...
          39.50 + 0.005*randn(n_cluster,1);
          16.4  + 0.4  *rand(n_sparse,1),  ...
          39.3  + 0.4  *rand(n_sparse,1)];
n_ps = size(lonlat, 1);

day = datenum(2021,9,14) + (0:12:12*186)';
n_ifg = numel(day);
years = (day - day(1)) / 365.25;

ph_mm = zeros(n_ps, n_ifg);
for k = 1:n_cluster
    ph_mm(k,:) = -15*years' + 0.4*randn(1, n_ifg);
end
for k = 1:n_sparse
    ph_mm(n_cluster+k,:) = (1 + 1.5*randn) * years' + 1.5*randn(1, n_ifg);
end

master_day = day(round(n_ifg/2));
lambda = 0.0555;
ref_ps = 1;
unwrap_ifg_index = 1:n_ifg;
ifg_list = 1:n_ifg;
bperp = randn(n_ifg, 1) * 50;

save(fullfile(tmp, 'ps_plot_ts_v-do.mat'), ...
     'ph_mm','lonlat','day','master_day','lambda', ...
     'ref_ps','unwrap_ifg_index','ifg_list','bperp','n_ps');

fprintf('Demo workdir : %s\n', tmp);
fprintf('Synthetic .mat ready. Launching PHASE_StaMPS...\n');

cd(fullfile(phase_root, 'PHASE_Preprocessing'));
app = PHASE_StaMPS;

pause(2);

app.PathEditField_2.Value = stamps_root;
app.PathEditField_3.Value = tmp;
app.MasterdateDatePicker.Value = datetime(2021,9,14);
app.ValueEditField.Value = 0.4;
app.TabGroup.SelectedTab = app.TSPointsTab;

fprintf('PHASE_StaMPS pre-configured. Now click the green button.\n');
