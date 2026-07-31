function applyConfig(engine, cfg)
%APPLYCONFIG Synchronize the modern form with the proven legacy engine.

isSEN = strcmp(cfg.constellation, 'SEN');
if isprop(engine,'GenerateCoherence'), engine.GenerateCoherence = cfg.generate_coherence; end
if isprop(engine,'GenerateLia'), engine.GenerateLia = cfg.generate_lia; end
engine.constellation = cfg.constellation;
engine.ConstellationSwitch.Value = ternary(isSEN, 'Sentinel1', 'CosmoSkyMed');
engine.Sentinel1Panel.Visible = onOff(isSEN);
engine.CosmoSkyMedPanel.Visible = onOff(~isSEN);
engine.StopFlag = false;

if isSEN
    engine.python_SEN = cfg.python;
    setPython(engine.PythonEnvironmentDropDown, ...
        engine.CustomPythonEnvironmentEditField, engine.PythonEnvironmentLabel, cfg.python);
    engine.update_processed_data_SEN = double(cfg.update_processed_data);
    engine.UpdateAlreadyProcessedDataCheckBox.Value = cfg.update_processed_data;
    engine.setSENUpdateControlsVisible(cfg.update_processed_data);
    engine.master_date_SEN = cfg.master_date;
    engine.MasterdateDatePicker.Value = asDate(cfg.master_date);
    engine.auto_master_SEN = double(cfg.auto_master);
    engine.AutoMasterCheckBox.Value = cfg.auto_master;
    engine.MasterdateDatePicker.Enable = onOff(~cfg.auto_master);
    engine.master_processing_SEN = double(~cfg.process_master);
    engine.MasterprocessingCheckBox.Value = cfg.process_master;
    engine.polarisation_SEN = cfg.polarisation;
    engine.PolarisationDropDown.Value = cfg.polarisation;
    engine.lon_min_SEN = cfg.lon_min; engine.MinlongitudeEditField.Value = cfg.lon_min;
    engine.lon_max_SEN = cfg.lon_max; engine.MaxlongitudeEditField.Value = cfg.lon_max;
    engine.lat_min_SEN = cfg.lat_min; engine.MinlatitudeEditField.Value = cfg.lat_min;
    engine.lat_max_SEN = cfg.lat_max; engine.MaxlatitudeEditField.Value = cfg.lat_max;
    % The beta performs all selected cleanup only after the engine returns.
    engine.slaves_removal_SEN = 1;
    engine.SlavesremovalafterprocessingCheckBox.Value = false;
    engine.dem_name_SEN = cfg.dem_name; engine.DEMinterferogramDropDown.Value = cfg.dem_name;
    engine.dem_file_SEN = cfg.dem_file; engine.DEMifgpathEditField.Value = cfg.dem_file;
    engine.dem_name_coreg_SEN = cfg.dem_name_coreg; engine.DEMcoregistrationDropDown.Value = cfg.dem_name_coreg;
    engine.dem_file_coreg_SEN = cfg.dem_file_coreg; engine.DEMcoregpathEditField.Value = cfg.dem_file_coreg;
    engine.dem_resampling_SEN = cfg.dem_resampling; engine.DEMresamplingmethodDropDown.Value = cfg.dem_resampling;
    engine.first_step_SEN = cfg.first_step; engine.FirststepDropDown.Value = num2str(cfg.first_step);
    generateTerrainProducts = cfg.generate_coherence || cfg.generate_lia;
    engine.coherence_tc_SEN = double(~generateTerrainProducts);
    engine.TerraincorrectedCoherenceandLIACheckBox.Value = generateTerrainProducts;
    engine.epsg_code_SEN = cfg.epsg_code; engine.EPSGcodeEditField.Value = cfg.epsg_code;
    engine.gptbin_path_SEN = cfg.gptbin_path; engine.PathEditField.Value = cfg.gptbin_path;
    engine.cpu_SEN = cfg.cpu; engine.CPUEditField.Value = cfg.cpu;
    engine.cache_SEN = cfg.cache; engine.CacheEditField.Value = cfg.cache;
    updateRoi(engine.roi_SEN, cfg);
else
    engine.python_CSK = cfg.python;
    setPython(engine.PythonEnvironmentDropDown_2, ...
        engine.CustomPythonEnvironmentEditField_2, engine.CustomPythonEnvironmentLabel, cfg.python);
    engine.master_date_CSK = cfg.master_date;
    engine.MasterdateDatePicker_2.Value = asDate(cfg.master_date);
    engine.auto_master_CSK = double(cfg.auto_master);
    engine.AutoMasterCheckBox_2.Value = cfg.auto_master;
    engine.MasterdateDatePicker_2.Enable = onOff(~cfg.auto_master);
    engine.master_processing_CSK = double(~cfg.process_master);
    engine.MasterprocessingCheckBox_2.Value = cfg.process_master;
    engine.lon_min_CSK = cfg.lon_min; engine.MinlongitudeEditField_2.Value = cfg.lon_min;
    engine.lon_max_CSK = cfg.lon_max; engine.MaxlongitudeEditField_2.Value = cfg.lon_max;
    engine.lat_min_CSK = cfg.lat_min; engine.MinlatitudeEditField_2.Value = cfg.lat_min;
    engine.lat_max_CSK = cfg.lat_max; engine.MaxlatitudeEditField_2.Value = cfg.lat_max;
    engine.slaves_removal_CSK = 1;
    engine.SlavesremovalafterprocessingCheckBox_2.Value = false;
    engine.dem_name_CSK = cfg.dem_name; engine.DEMinterferogramDropDown_2.Value = cfg.dem_name;
    engine.dem_file_CSK = cfg.dem_file; engine.DEMifgpathEditField_2.Value = cfg.dem_file;
    engine.num_gcp_CSK = cfg.num_gcp; engine.CoregistrationGCPsnumberEditField.Value = cfg.num_gcp;
    engine.first_step_CSK = cfg.first_step; engine.FirststepDropDown_2.Value = num2str(cfg.first_step);
    generateTerrainProducts = cfg.generate_coherence || cfg.generate_lia;
    engine.coherence_tc_CSK = double(~generateTerrainProducts);
    engine.TerraincorrectedCoherenceandLIACheckBox_2.Value = generateTerrainProducts;
    engine.epsg_code_CSK = cfg.epsg_code; engine.EPSGcodeEditField_2.Value = cfg.epsg_code;
    engine.gptbin_path_CSK = cfg.gptbin_path; engine.PathEditField_2.Value = cfg.gptbin_path;
    engine.cpu_CSK = cfg.cpu; engine.CPUEditField_2.Value = cfg.cpu;
    engine.cache_CSK = cfg.cache; engine.CacheEditField_2.Value = cfg.cache;
    updateRoi(engine.roi_CSK, cfg);
end
end

function setPython(dropdown, customField, customLabel, value)
standard = {'python','python3','python3.11'};
if any(strcmp(value, standard)) && any(strcmp(value, dropdown.Items))
    dropdown.Value = value;
    customLabel.Visible = 'off'; customField.Visible = 'off';
else
    dropdown.Value = 'Other';
    customField.Value = value;
    customLabel.Visible = 'on'; customField.Visible = 'on';
end
end

function value = asDate(textValue)
try
    value = datetime(textValue, 'InputFormat', 'yyyyMMdd');
catch
    value = datetime('today');
end
end

function updateRoi(roi, cfg)
try
    if ~isempty(roi) && isvalid(roi)
        roi.Position = [cfg.lat_min, cfg.lon_min, ...
            cfg.lat_max - cfg.lat_min, cfg.lon_max - cfg.lon_min];
    end
catch
end
end

function value = onOff(condition)
if condition, value = 'on'; else, value = 'off'; end
end

function value = ternary(condition, yesValue, noValue)
if condition, value = yesValue; else, value = noValue; end
end
