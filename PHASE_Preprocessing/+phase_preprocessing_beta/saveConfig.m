function pathValue = saveConfig(rootDir, cfg)
%SAVECONFIG Write the same input MAT contract used by the stable MLAPP.

folder = fullfile(rootDir, 'PHASE_Preprocessing');
if ~isfolder(folder), mkdir(folder); end
pathValue = fullfile(folder, 'input_preprocessing.mat');

constellation = cfg.constellation;
python = cfg.python;
master_date = cfg.master_date;
auto_master = double(cfg.auto_master);
master_processing = double(~cfg.process_master);
lon_min = cfg.lon_min; lat_min = cfg.lat_min;
lon_max = cfg.lon_max; lat_max = cfg.lat_max;
slaves_removal = double(~cfg.remove_slaves_after_processing);
remove_split_after_processing = double(cfg.remove_split_after_processing);
remove_coreg_after_processing = double(cfg.remove_coreg_after_processing);
remove_ifg_after_processing = double(cfg.remove_ifg_after_processing);
dem_name = cfg.dem_name; dem_file = cfg.dem_file;
first_step = cfg.first_step;
auto_epsg = double(cfg.auto_epsg);
generate_coherence = double(cfg.generate_coherence);
generate_lia = double(cfg.generate_lia);
coherence_tc = double(~(cfg.generate_coherence || cfg.generate_lia));
epsg_code = cfg.epsg_code; gptbin_path = cfg.gptbin_path;
cpu = cfg.cpu; cache = cfg.cache;

if strcmp(cfg.constellation, 'SEN')
    polarisation = cfg.polarisation;
    dem_name_coreg = cfg.dem_name_coreg;
    dem_file_coreg = cfg.dem_file_coreg;
    dem_resampling = cfg.dem_resampling;
    save(pathValue, 'constellation', 'python', 'master_date', 'auto_master', ...
        'master_processing', 'polarisation', 'lon_min', 'lat_min', 'lon_max', ...
        'lat_max', 'slaves_removal', 'dem_name', 'dem_file', 'dem_name_coreg', ...
        'dem_file_coreg', 'dem_resampling', 'first_step', 'coherence_tc', ...
        'auto_epsg', 'generate_coherence', 'generate_lia', ...
        'remove_split_after_processing', 'remove_coreg_after_processing', ...
        'remove_ifg_after_processing', ...
        'epsg_code', 'gptbin_path', 'cpu', 'cache', '-mat');
else
    num_gcp = cfg.num_gcp;
    save(pathValue, 'constellation', 'python', 'master_date', 'auto_master', ...
        'master_processing', 'lon_min', 'lat_min', 'lon_max', 'lat_max', ...
        'slaves_removal', 'dem_name', 'dem_file', 'first_step', 'num_gcp', ...
        'coherence_tc', 'auto_epsg', 'generate_coherence', 'generate_lia', ...
        'remove_split_after_processing', 'remove_coreg_after_processing', ...
        'remove_ifg_after_processing', ...
        'epsg_code', 'gptbin_path', 'cpu', 'cache', '-mat');
end
end
