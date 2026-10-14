function applyConfig(engine, cfg)
%APPLYCONFIG Set the headless preprocessing runner's execution options.
engine.constellation = string(cfg.constellation);
engine.StopFlag = false;
engine.GenerateCoherence = logical(cfg.generate_coherence);
engine.GenerateLia = logical(cfg.generate_lia);
engine.update_processed_data_SEN = double(cfg.update_processed_data);
end
