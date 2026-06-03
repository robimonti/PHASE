% launch_ts_picker_demo.m
%
% Genera un ps_plot_ts_v-do.mat sintetico in una tempdir e lancia
% ts_export_picker. Serve solo per ispezionare la GUI senza un
% workdir PSI reale.

% Derive paths from this script's location so the demo runs on any
% checkout, not just F:/phase.
this_dir   = fileparts(mfilename('fullpath'));
phase_root = fileparts(this_dir);
addpath(fullfile(phase_root, 'StaMPS', 'matlab'));
addpath(fullfile(phase_root, 'StaMPS', 'matlab_compat'));

tmp = tempname;
mkdir(tmp);

% PS sintetici sull'AOI Calabria: 600 PS in cluster + sparsi
rng(42);
n_cluster = 200;
n_sparse  = 400;
lonlat_cluster = [16.60 + 0.005*randn(n_cluster,1), ...
                  39.50 + 0.005*randn(n_cluster,1)];
lonlat_sparse  = [16.4 + 0.4*rand(n_sparse,1), ...
                  39.3 + 0.4*rand(n_sparse,1)];
lonlat = [lonlat_cluster; lonlat_sparse];
n_ps   = size(lonlat,1);

day = datenum(2021,9,14) + (0:12:12*186)';   % 187 acquisitions, 12-day cadence
n_ifg = numel(day);

% Cluster mostra subsidence -15 mm/yr; sparse rumore ±2 mm/yr
years = (day - day(1)) / 365.25;
ph_mm = zeros(n_ps, n_ifg);
for k = 1:n_cluster
    ph_mm(k,:) = -15*years' + 0.4*randn(1,n_ifg);
end
for k = 1:n_sparse
    ph_mm(n_cluster+k,:) = (1 + 1.5*randn) * years' + 1.5*randn(1,n_ifg);
end

master_day = day(round(n_ifg/2));
lambda = 0.0555;
ref_ps = 1;
unwrap_ifg_index = 1:n_ifg;
ifg_list = 1:n_ifg;
bperp = randn(n_ifg,1)*50;

matfile = fullfile(tmp, 'ps_plot_ts_v-do.mat');
save(matfile, 'ph_mm','lonlat','day','master_day','lambda', ...
              'ref_ps','unwrap_ifg_index','ifg_list','bperp','n_ps');

fprintf('Synthetic .mat written: %s\n', matfile);
fprintf('Launching ts_export_picker...\n');
ts_export_picker(tmp);
