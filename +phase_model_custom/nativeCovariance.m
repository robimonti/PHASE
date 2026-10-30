function settings=nativeCovariance(residuals_grid_full,xyIN_AOI,t_relIN,config)
% Native empirical product covariance and native analytical model fits.
addpath(fullfile(phase_model_beta.projectRoot(),'MatlabFunctions'));
detrend_method='cleanObs';
dtCov_STC2D=config.(['dtCov_STC' config.projDim]);
dsCov_STC2D=config.(['dsCov_STC' config.projDim]);
tCovModel_STC2D=config.(['tCovModel_STC' config.projDim]);
sCovModel_STC2D=config.(['sCovModel_STC' config.projDim]);
if isnan(dsCov_STC2D)
    distances=pdist(xyIN_AOI); distances=distances(distances>0);
    dsCov_STC2D=median(distances)/3;
end
if ~isfinite(dsCov_STC2D)||dsCov_STC2D<=0, error('PHASE:customCovarianceBins','Spatial covariance needs distinct locations and positive bin size.'); end
    % 3.7) Global Empirical Covariance Estimation
    disp('Computing global empirical covariance...');
    
        % 3.7.1) Covariance in TIME
        if isnan(dtCov_STC2D)
            dtCov = mode(diff(t_relIN)); if dtCov==0, dtCov=1; end
        else
            dtCov = dtCov_STC2D;
        end
        max_lag_time = ceil((t_relIN(end)-t_relIN(1))/dtCov) + 1;
        
        num_ps = size(residuals_grid_full, 1);
        accum_T_cells = cell(num_ps, 1);
        var_T_cells   = cell(num_ps, 1);   % store variance parts
        
        for i = 1:num_ps
            ts = residuals_grid_full(i, :);
            valid = ~isnan(ts);
            if sum(valid) > 1
                % Compute emp cov
                [local_acc, s_sq, n_p] = f1DEmpCovEst_accum(ts(valid), t_relIN(valid), dtCov, false);
                accum_T_cells{i} = local_acc;
                var_T_cells{i}   = [s_sq, n_p];
            end
        end
        
        % 3.7.1a) Aggregate variance (lag 0)
        Var_Global_T = sum(cat(1, var_T_cells{:}), 1);   % [SumSq, TotalCount]
        if numel(Var_Global_T) < 2 || Var_Global_T(2) == 0
            Variance_T = 1e-3; % Safe fallback
        else
            Variance_T = Var_Global_T(1) / Var_Global_T(2);
        end
        
        % 3.7.1b) Aggregate covariance (lag > 0)
        Global_T = zeros(max_lag_time, 3);
        for i = 1:num_ps
            acc = accum_T_cells{i};
            if ~isempty(acc)
                n = size(acc, 1);
                if n > size(Global_T, 1), Global_T(n, 3) = 0; end
                Global_T(1:n, :) = Global_T(1:n, :) + acc;
            end
        end
        
        valid_t = Global_T(:, 2) > 0;
        lag_t_raw    = Global_T(valid_t, 3) ./ Global_T(valid_t, 2);
        eCovF_mean_t = Global_T(valid_t, 1) ./ Global_T(valid_t, 2);
        Cecf_mean_t  = 1 ./ Global_T(valid_t, 2);
        
        % 3.7.1c) Combine: prepend variance at lag 0
        lag_t        = [0; lag_t_raw];
        eCovF_mean_t = [Variance_T; eCovF_mean_t];
        Cecf_mean_t  = [1/Var_Global_T(2); Cecf_mean_t];   % Weight for variance
        
        % Sort
        [lag_t, idx] = sort(lag_t);
        eCovF_mean_t = eCovF_mean_t(idx);
        Cecf_mean_t  = Cecf_mean_t(idx);
    
    
        % 3.7.2) Covariance in SPACE
        if isnan(dsCov_STC2D)
            dsCov = prctile(sqrt((xyIN_AOI(:,1)).^2 + (xyIN_AOI(:,2)).^2), 25);
        else
            dsCov = dsCov_STC2D;
        end 
        max_dist = sqrt(range(xyIN_AOI(:,1))^2 + range(xyIN_AOI(:,2))^2);
        edges = 0:dsCov:ceil(max_dist / dsCov) * dsCov;
        n_epochs = size(residuals_grid_full, 2);
        accum_S_cells = cell(n_epochs, 1);
        var_S_cells   = cell(n_epochs, 1);
        var_S_cells1   = cell(n_epochs, 1);
    
        % Broadcast variable for parfor
        method_flag = detrend_method; 
        
        for i = 1:n_epochs
            sp = residuals_grid_full(:, i);
            valid = ~isnan(sp);
            
            if sum(valid) > 1
                sp_val = sp(valid);
                coords_val = xyIN_AOI(valid, :);
                
                % Conditional centering
                switch method_flag
                    case 'cleanObs'
                    % In cleanObs mode, the global surface removal might leave 
                    % epoch-wise offsets. The mean is removed locally strictly 
                    % for covariance estimation to ensure the tail decays to zero.
                    sp_input = sp_val - mean(sp_val);
    
                    [~, s_sq1, n_p1] = f2DEmpCovEst_accum(sp_val, coords_val, edges);
                    var_S_cells1{i}   = [s_sq1, n_p1];
                    case 'residualAtm'
                    % In residualAtm mode, the iterative LSC has already removed 
                    % the atmosphere (mean/ramp).
                    sp_input = sp_val;
                end
                
                [local_acc, s_sq, n_p] = f2DEmpCovEst_accum(sp_input, coords_val, edges);
                accum_S_cells{i} = local_acc;
                var_S_cells{i}   = [s_sq, n_p];
            end
        end
        
        % 3.7.2a) Aggregate variance (lag 0)
        switch method_flag
            case 'cleanObs'
                Var_Global_S = sum(cat(1, var_S_cells1{:}), 1);
            case 'residualAtm'
                Var_Global_S = sum(cat(1, var_S_cells{:}), 1);
        end
        if numel(Var_Global_S) < 2 || Var_Global_S(2) == 0
            Variance_S = 1e-3; % Safe fallback
        else
            Variance_S = Var_Global_S(1) / Var_Global_S(2);
        end
        
        % 3.7.2b) Aggregate covariance (lag > 0)
        Global_S = sum(cat(3, accum_S_cells{:}), 3);
        
        valid_s = Global_S(:, 2) > 0;
        lag_s_raw    = Global_S(valid_s, 3) ./ Global_S(valid_s, 2);
        eCovF_mean_s = Global_S(valid_s, 1) ./ Global_S(valid_s, 2);
        Cecf_mean_s  = 1 ./ Global_S(valid_s, 2);
        
        % 3.7.2c) Combine
        lag_s        = [0; lag_s_raw];
        eCovF_mean_s = [Variance_S; eCovF_mean_s];
        Cecf_mean_s  = [1/Var_Global_S(2); Cecf_mean_s];
        
        [lag_s, idx] = sort(lag_s);
        eCovF_mean_s = eCovF_mean_s(idx);
        Cecf_mean_s  = Cecf_mean_s(idx);
        
        disp('Global covariance estimation completed.');

    
    % 3.8) Fit TIME covariance
    x_data_t = lag_t;
    if length(lag_t) < 20
        lag_t1 = linspace(0,lag_t(end),101);
    else
        lag_t1 = lag_t;
    end
    y_data_t = eCovF_mean_t;
    w_data_t = Cecf_mean_t;
    
    % Use the chosen model
    switch tCovModel_STC2D
        case 'exponential'
            [c_opt_t, mCovF1] = fit_exp_auto(x_data_t, y_data_t, w_data_t);
        case 'gaussian'
            [c_opt_t, mCovF1] = fit_gaussian_auto(x_data_t, y_data_t, w_data_t);
        case {'gaussian_cos', 'gaussian_cosine', 'gaussian+cosine'}
            [c_opt_t, mCovF1] = fit_gaussian_cosine_auto(x_data_t, y_data_t, w_data_t);
    end
    c1_t = c_opt_t;
    
    
    % 3.9) Fit SPACE covariance
    x_data_s = lag_s;
    if length(lag_s) < 20
        lag_s1 = linspace(0,lag_s(end),101);
    else
        lag_s1 = lag_s;
    end
    y_data_s = eCovF_mean_s;
    w_data_s = Cecf_mean_s;
    
    % Use the chosen model
    switch sCovModel_STC2D
        case 'exponential'
            [c_opt_s, mCovF2] = fit_exp_auto(x_data_s, y_data_s, w_data_s);
        case 'gaussian'
            [c_opt_s, mCovF2] = fit_gaussian_auto(x_data_s, y_data_s, w_data_s);
        case {'gaussian_cos', 'gaussian_cosine', 'gaussian+cosine'}
            [c_opt_s, mCovF2] = fit_gaussian_cosine_auto(x_data_s, y_data_s, w_data_s);
    end
    c1_s = c_opt_s;
    
    
    % 3.10) Calculate noise (nugget)
    % Noise is the difference between Empirical(0) and Model(0)
        c1_m = mean([c1_t(1), c1_s(1)]);
        c1_t(1) = c1_m; c1_s(1) = c1_m;
    C_signal_at_0 = mCovF2(c1_s, 0);
    Empirical_at_0 = mean([eCovF_mean_t(1), eCovF_mean_s(1)]);
    
    A_noise = Empirical_at_0 - C_signal_at_0;
    
        % Safety: noise cannot be negative. 
        % If signal model is higher than data, assume zero noise.
        if A_noise < 0
            A_noise = 1e-5; 
        end


settings=struct('dimension',config.projDim,'native',true,'timeParameters',c1_t, ...
    'spaceParameters',c1_s,'timeFunction',mCovF1,'spaceFunction',mCovF2, ...
    'noiseVariance',A_noise,'signalVariance',c1_m,'lagTime',lag_t, ...
    'covarianceTime',eCovF_mean_t,'lagSpace',lag_s,'covarianceSpace',eCovF_mean_s, ...
    'centeringPolicy','Original cleanObs: spatial epoch centering at nonzero lags; raw second moment at lag zero; temporal products unchanged');
end
