function cfg = defaultConfig()
%DEFAULTCONFIG Stable PHASE preprocessing defaults in explicit form.

cfg = struct();
cfg.constellation = 'SEN';
cfg.python = 'python';
cfg.update_processed_data = false;
cfg.master_date = '20200722';
cfg.auto_master = true;
cfg.process_master = true;
cfg.polarisation = 'VV';
cfg.lon_min = -180;
cfg.lat_min = -90;
cfg.lon_max = 180;
cfg.lat_max = 90;
cfg.remove_slaves_after_processing = false;
cfg.remove_split_after_processing = false;
cfg.remove_coreg_after_processing = false;
cfg.remove_ifg_after_processing = false;
cfg.dem_name = 'SRTM 1Sec HGT';
cfg.dem_file = '';
cfg.dem_name_coreg = 'SRTM 1Sec HGT';
cfg.dem_file_coreg = '';
cfg.dem_resampling = 'NEAREST_NEIGHBOUR';
cfg.first_step = 1;
cfg.auto_epsg = true;
cfg.epsg_code = 32631;
cfg.generate_coherence = true;
cfg.generate_lia = true;
cfg.gptbin_path = 'C:\Program Files\snap\bin\gpt';
cfg.cpu = 8;
cfg.cache = '26G';
cfg.num_gcp = 10000;
end
