function [lonlatIN_AOI_DET1D, dates_full, t_full, final_signal_orig, final_signal_out, final_std_out] = ...
    STmodel_DET1D(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, centerline_data, figsDir, minMonths, ...
    gS_input_path, gS_output_path, gS_job_path, gS_synth_path, detectedOS, varNoise_method, varNoise_manual, num_spl_method, ...
    spline_method, num_spl_row_manual, num_spl_col_manual, lambda_method, lambda_manual, utmZone, xyIN_AOI, step_t, ...
    detrend_method, use_inclined_means, poly_degree, markerSize_DET1D)

% STmodel_DET1D Performs spatio-temporal deterministic modelling for displacement data
%
% Inputs:
%   displIN_AOI           - n x m matrix of displacement data (n points, m epochs)
%   PSidIN_AOI            - n x 1 vector of point IDs
%   t_dateIN              - m x 1 vector of absolute dates (datetime)
%   t_relIN               - m x 1 vector of relative times (days)
%   centerline_data       - Structure with fields:
%                          - xy_centerline: k x 2 matrix of centerline coordinates
%                          - s: k x 1 vector of cumulative distances along centerline
%                          - projected_distances: n x 1 vector of projected distances
%   figsDir               - Directory path for saving figures
%   file_out              - Minimum period length for splines interpolation
%   gS_input_path         - path to input folder for geoSplinter
%   gS_output_path        - path to output folder for geoSplinter
%   gS_job_path           - path to job folder for geoSplinter
%   gS_synth_path         - path to synthesis folder for geoSplinter
%   detectedOS            - detected operative system
%   varNoise              - a-priori noise variance method 'auto'/'manual'
%   varNoise_manual       - a-priori noise variance value for manual selection (otherwise NaN)
%   num_spl_method        - method for splines number selection 'auto'/'manual'
%   spline_method         - method for 'auto' splines selection 'variance'/'MDL'/'F-test'/'chi2_test'
%   num_spl_row_manual    - manual number of splines for rows (otherwise NaN)
%   num_spl_col_manual    - manual number of splines for columns (otherwise NaN)
%   lambda_method         - method for lambda estimation 'auto'/'manual'
%   lambda_manual         - value of lambda for 'manual' case' (otherwise NaN)
%   utmZone               - UTM zone for reprojection of coordinates
%   xyIN_AOI              - coordinates of PS inside the AOI
%   step_t                - time step for the estimation
%   detrend_method        - option to select the processing method: residual atmospheric errors (residualAtm) or clean data (cleanObs)
%   use_inclined_means    - flag for using flat or inclines plane for residual atmospheric removal
%   poly_degree           - maximum polynomial degree
%   markerSize_STC2D      - marker size for the GIF plot
% 
%
% Output:
%   lonlatIN_AOI_DET1D    - lon/lat coordinates of estimated centerline points
%   dates_full            - dates prediction
%   t_full                - relative time prediction
%   final_signal_orig     - modelled signal in observation coordinates
%   final_signal_out      - modelled signal in query coordinates
%   final_std_out         - uncertainty

% validate centerline inputs
if isempty(centerline_data) || isempty(centerline_data.xy_centerline) || isempty(centerline_data.s) || isempty(centerline_data.projected_distances)
    error('centerline_data or its fields are empty for 1D interpolation.');
end


disp('--- 1D Deterministic-based modelling started... ---');


% interpolate along precomputed centerline
    % grid_x = centerline_data.xy_centerline(:,1);
    % grid_y = centerline_data.xy_centerline(:,2);
s_full_SAR = centerline_data.s;
s_SAR = centerline_data.projected_distances;

% sorting index based on position along centerline
[~, idx_SAR] = sort(s_SAR);




% identify the master epoch (date where displacements are zero)
tolerance = 1e-6;    % small tolerance for numerical precision
master_idx = find(all(abs(displIN_AOI) < tolerance, 1));
if isempty(master_idx)
    warning('No master epoch found where displacements are nearly zero. Using first epoch.');
    master_idx = 1;
elseif length(master_idx) > 1
    % select epoch closest to mean of t_relIN
    [~, closest_idx] = min(abs(t_relIN(master_idx) - mean(t_relIN)));
    master_idx = master_idx(closest_idx);
    warning('Multiple near-zero displacement epochs found. Using epoch at t_relIN = %.2f.', t_relIN(master_idx));
end
master_t_rel = t_relIN(master_idx);




%% 1) SPLINES-BASED OUTLIER REJECTION
disp('======== Step 1 ========');
disp('Splines-based outlier rejection started...');

% minimum time interval in days
minT = minMonths * 30;  % months -> days

% define basis functions from geoSplinter
phi_3 = @(csi) (csi >= 0 & csi < 1) .* ((2 - csi).^3 - 4 * (1 - csi).^3) / 6 + ...
               (csi >= 1 & csi < 2) .* ((2 - csi).^3) / 6;
phi_4 = @(csi) (csi >= 0 & csi < 2) .* ((2 - csi).^3) / 6;

% initialize result storage variables
data_reg_spl = zeros(size(displIN_AOI'));
t_reg_spl = zeros(size(displIN_AOI));
data_raw_cleanT = zeros(size(displIN_AOI));

% define folder for geoSplinter based on current OS
switch detectedOS
    case 'windows'
        % Windows case
        gS_dir = (strcat('.', filesep, 'geoSplinter', filesep, 'windows')); 

    case 'linux'
        % Linux case
        gS_dir = (strcat('.', filesep, 'geoSplinter', filesep, 'linux')); 

    case 'macos_intel'
        % macOS (Intel) case
        gS_dir = (strcat('.', filesep, 'geoSplinter', filesep, 'macos_intel')); 

    case 'macos_apple_silicon'
        % macOS (Apple Silicon) case
        gS_dir = (strcat('.', filesep, 'geoSplinter', filesep, 'macos_apple_silicon')); 

    otherwise
        error('Unsupported operating system.');
end

% process each PS individually
for id = 1:size(displIN_AOI,1)
    
        fprintf('Currently processing PS %i\n', PSidIN_AOI(id));  % display current PS ID
        
        % set data and time variables for current PS
        data_reg = displIN_AOI';                      % full dataset, transposed
        t_reg = t_relIN';                              % time vector
        data_reg_PS = data_reg(:,id);                 % data for current PS
        
        % write input data to file for geoSplinter
        gS_filename = strcat('PS_', num2str(PSidIN_AOI(id)));
        writematrix([t_reg, data_reg_PS], fullfile(gS_input_path, [gS_filename, '.txt']), 'Delimiter', 'space');
        ext = 'cub';                                  % extension for cubic splines

        % define job file parameters
        data_dim = 1;                                 % data dimension
        type_spl = 2;                                 % spline type (1 = linear, 2 = cubic)
        file_inp = strcat(gS_filename, '.txt');       % input filename
        file_out = strcat(gS_filename, '_', ext);     % output filename
        num_obs = size(data_reg_PS,1);                % number of observations
        t_in = t_reg(1);                              % start time
        t_fin = round(t_reg(end));                    % end time
        lambda = 0;                                   % regularization parameter
        num_sig = 10;                                 % number of significant digits
        
        % define MDL-based spline search parameters
        max_n_spl = round(t_reg(end) / minT, 0);
        min_n_spl = 3;
        num_spl_candidates = max_n_spl - min_n_spl + 1; 

        % initialize variables
        data_spline_tmp = NaN(size(displIN_AOI, 2), num_spl_candidates);
        MDL_tmp = NaN(num_spl_candidates, 1);
    
        % initialize MDL optimization
        k = 1;                                        % iteration counter for MDL search
        for i = min_n_spl:max_n_spl
    
            num_spl = i;                              % candidate number of splines

            % write the job file
            jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, num_spl, ...
                    t_in, t_fin, lambda, num_sig, gS_input_path, gS_output_path, gS_job_path);
    
            % run geoSplinter_analysis with the job file
            job_execution = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_analysis'), fullfile('.', gS_job_path, [file_out, '.job']));
            status = system(job_execution);
            if status ~= 0
                error('Error executing geoSplinter_analysis for file: %s', file_out);
            end

            % define the output file base name from geoSplinter run
            mat_file = fullfile('.', gS_output_path, strcat(file_out, '.mat.txt'));
            hdr_file = fullfile('.', gS_output_path, strcat(file_out, '.hdr.txt'));

            % import normal matrix and normal vector
            fid = fopen(mat_file, 'r');
            if fid == -1
                error('Cannot open file: %s', mat_file);
            end

            n_params = num_spl + 2; % total parameters (splines + 2)
            normal_matrix = zeros(n_params, n_params);
            normal_vector = zeros(n_params, 1);
            
            while ~feof(fid)
                line = fgetl(fid);
                if ~ischar(line) || isempty(strtrim(line))
                    continue;
                end
                
                % extract Normal Matrix (not "as stored")
                if contains(line, 'Normal Matrix') && ~contains(line, 'as stored') && ~contains(line, 'Diagonal')
                    r = 1;
                    while r <= n_params
                        line = fgetl(fid);
                        if ~ischar(line) || isempty(strtrim(line))
                            continue;
                        end
                        % parse decimal format
                        values = sscanf(line, repmat('%f ', 1, n_params));
                        if length(values) == n_params
                            normal_matrix(r, :) = values';
                            r = r + 1;
                        else
                            error('Unexpected format in Normal Matrix row %d: expected %d values, got %d', ...
                                  r, n_params, length(values));
                        end
                    end
                end
                
                % extract Normal Vector
                if contains(line, 'Normal Vector')
                    r = 1;
                    while r <= n_params
                        line = fgetl(fid);
                        if ~ischar(line) || isempty(strtrim(line))
                            continue;
                        end
                        normal_vector(r) = str2double(line);
                        r = r + 1;
                    end
                end

                if contains(line, 'Sigma Zero')
                    tokens = regexp(line, 'Sigma Zero\s*:\s*([\d.]+)', 'tokens');
                    if ~isempty(tokens)
                        sigma0_est = str2double(tokens{1}{1});
                    end
                    break;
                end
            end
            fclose(fid);
            
            % import knots
            fid = fopen(hdr_file, 'r');
            if fid == -1
                error('Cannot open file: %s', hdr_file);
            end
            
            xMin = []; xMax = []; deltaX = [];
            while ~feof(fid)
                line = fgetl(fid);
                if ~ischar(line) || isempty(strtrim(line))
                    continue;
                end
                if contains(line, 'xMin')
                    xMin = sscanf(line, 'xMin   :   %e');
                elseif contains(line, 'xMax')
                    xMax = sscanf(line, 'xMax   :   %e');
                elseif contains(line, 'deltaX')
                    deltaX = sscanf(line, 'deltaX :   %e');
                end
            end
            fclose(fid);
            
            % construct knots
            knots = xMin:deltaX:xMax;

            % apply regularization for master epoch constraint
            % build design matrix A
            A = zeros(num_obs, n_params);
            for n = 1:num_obs
                t = t_reg(n);
                i_i = find(t >= knots(1:end-1) & t < knots(2:end), 1) - 1;
                if isempty(i_i)
                    if t <= knots(1)
                        i_i = 0;
                    elseif t >= knots(end)
                        i_i = length(knots) - 2;
                    end
                end
                csi_i = (t - knots(i_i + 1)) / deltaX;
                
                alpha = zeros(4, 1);
                alpha(1) = phi_4(csi_i + 1);  % B_{i_i-1}
                alpha(2) = phi_3(csi_i);      % B_{i_i}
                alpha(3) = phi_3(1 - csi_i);  % B_{i_i+1}
                alpha(4) = phi_4(2 - csi_i);  % B_{i_i+2}
                
                for s = -1:2
                    j = i_i + s + 1;
                    if j >= 1 && j <= n_params
                        A(n, j) = alpha(s + 2);
                    end
                end
            end
            
            % set B from A
            B = A(master_idx, :);
            
            % regularization
            R = B' * B;
            lambdaR = 10;
            N_modified = normal_matrix + lambdaR * R;
            x_constrained = N_modified \ normal_vector;
            
            % compute estimated observations
            estimated_obs = A * x_constrained;
            clear normal_matrix normal_vector;
            
            % calculate MDL for current spline number
            SSR = sum((data_reg_PS - estimated_obs).^2);           % sum of squared residuals
            data_spline_tmp(:,k) = estimated_obs;                  % temporary storage for current spline solution
            MDL_tmp(k) = log(num_obs) * num_spl / num_obs + log(SSR);  % MDL calculation
    
            k = k + 1;  % update iteration counter
    
        end
    
        % identify best num_spl that minimizes MDL index
        [~, MDL_idx] = min(MDL_tmp);
        
        % store final spline results and residuals
        data_spline = data_spline_tmp(:,MDL_idx);             % best spline fit
        res_det = data_reg_PS - data_spline;                  % residuals for the best spline n°
        clear data_spline_tmp MDL_tmp;

        Q = eye(size(res_det,2),size(res_det,2));
        N = A' * (Q \ A);

        % outlier removal
        Z_lim = norminv(1-0.05/2);
        Z_obs = abs(res_det ./ sqrt(diag(sigma0_est * (Q - A / N * A'))));
        clear Q N;
        data_reg_PS(Z_obs > Z_lim) = NaN;
        t_reg(Z_obs > Z_lim) = NaN;

        % store cleaned data and time vector for this PS
        data_raw_cleanT(id,:) = data_reg_PS';
        t_reg_spl(id,:) = t_reg';
        data_reg_spl(:,id) = data_spline;
        
        fclose('all');
        pause(0.01);
end

disp('Splines-based outlier rejection completed.');




%% 2) PREPARE VARIABLES
disp('======== Step 2 ========');
disp('Variables preparation started...');

% Define variables
displIN_AOI_sorted = displIN_AOI(idx_SAR, :);
data_reg_spl = data_reg_spl';
data_reg_spl_sorted = data_reg_spl(idx_SAR, :);
dataRawCleanT_sorted = data_raw_cleanT(idx_SAR, :);
t_reg_sorted = t_reg_spl(idx_SAR, :);

% Compute residuals of observations - splines interpolation
res_spl = displIN_AOI_sorted - data_reg_spl_sorted;

% Reorder s_SAR positions along centerline based on increasing order
s_SAR_orig = s_SAR;
s_SAR = s_SAR(idx_SAR);

% Define dimensions
np = length(s_SAR);
nt = length(t_relIN);

% Reshape cleaned observations to column vector (ordered by time, then point)
dataRawCleanT_sorted_col = reshape(dataRawCleanT_sorted', [], 1);
dataRawCleanT_sorted_col(isnan(dataRawCleanT_sorted_col)) = [];

% Create the full grid of observed positions (s) and times
grid_t0 = repmat(t_relIN, np, 1);  % [np x nt]
grid_s0 = repmat(s_SAR, 1, nt);    % [np x nt]

% Filter out missing observations
mask_obs = isnan(dataRawCleanT_sorted);  % [np x nt]
grid_t_flt = grid_t0;
grid_t_flt(mask_obs) = NaN;
grid_s_flt = grid_s0;
grid_s_flt(mask_obs) = NaN;

% Make them column vectors (same ordering: time-major, then point)
mask_obs_col = reshape(mask_obs', [], 1);  % [nt * np x 1]
grid_t_flt_col = reshape(grid_t_flt', [], 1);
grid_t_flt_col(isnan(grid_t_flt_col)) = [];
grid_s_flt_col = reshape(grid_s_flt', [], 1);
grid_s_flt_col(isnan(grid_s_flt_col)) = [];

% Do the same for the full grid without removed nans
t_full_0 = reshape(grid_t0', [], 1);  % [nt * np x 1]
s_full_0 = reshape(grid_s0', [], 1);

residuals_grid_full = dataRawCleanT_sorted;

disp('Variables preparation completed.');




%% 3) LEAST SQUARES MODEL REDUCTION
disp('======== Step 3 ========');
disp('Least Squares model reduction started...');

% observations
y0_lsc = reshape(residuals_grid_full', [], 1);
y0_lsc(isnan(y0_lsc)) = [];

% Prepare the environment
% centering of coordinates to fix ill-conditioning
s0 = mean(grid_s_flt_col);
t0 = mean(grid_t_flt_col);
    
    % centered versions for fitting (observed points)
    s_centered = grid_s_flt_col - s0;
    t_centered = grid_t_flt_col - t0;
    % centered versions for reconstruction on the full regular grid
    s_full_centered = s_full_0 - s0;
    t_full_centered = t_full_0 - t0;


% 3.1) Covariance matrix
alpha_sig = 0.05;
variances = ones(length(y0_lsc), 1);
idx_master = (grid_t_flt_col == master_t_rel);
if any(idx_master)
    variances(idx_master) = 0.1;
end
Q_ll = spdiags(variances, 0, length(y0_lsc), length(y0_lsc));


% 3.2) Design matrix construction
n_epochs = length(t_relIN);

switch detrend_method

    case 'residualAtm'
    % 3.2.1) Residual errors in the data (Orbits + Atmosphere)

    % Global orbital trend (static plane a*s)
    A_trend = s_centered;
    n_trend_cols = size(A_trend, 2);
    
    if ~use_inclined_means
        % 3.2.2a Small area: constant atmosphere per epoch
       
        n_PS = size(residuals_grid_full, 1);
        
        % Identity matrix for epochs, replicated for each PS
        A_atmos_full = repmat(eye(n_epochs), n_PS, 1);
        
        % Extract valid observations
        A_atmos = A_atmos_full(~mask_obs_col, :);
        
    else
        % 3.2.2b) Large area: ramped atmosphere per epoch
        
        % 2 parameters per epoch: [Bias, SlopeS]
        % Total atmos params = 2 * n_epochs
        
        % For each observation (i, t) -> columns corresponding to epoch t.
        % The row is: [ ... 0 0 1 s_i 0 0 ... ] at the columns for epoch t.
        
        % Map each valid observation to its epoch index (1 to n_epochs)
        [~, obs_epoch_idx] = ismember(grid_t_flt_col, t_relIN);
        
        % Base indices for the 2 parameters of each epoch
        % Cols for epoch k are: (k-1)*2 + [1, 2]
        col_bias  = (obs_epoch_idx - 1) * 2 + 1;
        col_slopeS = (obs_epoch_idx - 1) * 2 + 2;
        
        n_total_obs = length(grid_t_flt_col);
        n_atmos_params = 2 * n_epochs;
        
        % Create sparse matrix (memory efficient for large N)
        % Rows: 1 to N, Cols: computed above, Values: 1, x
        I = repmat((1:n_total_obs)', 2, 1);
        J = [col_bias; col_slopeS];
        V = [ones(n_total_obs, 1); s_centered];
        
        A_atmos = sparse(I, J, V, n_total_obs, n_atmos_params);
        
        % Convert to full if N is small enough to speed up solvers
        if n_total_obs < 10000
            A_atmos = full(A_atmos);
        end

    end
    
    % Combine Global Trend + Atmosphere
    A_full = [A_trend, A_atmos];


    % 3.3) LS with backward elimination
    % Track which parameters are currently kept in the model (initially, everything is kept)
    mask_params = true(size(A_full, 2), 1);
    
    % Loop until no more parameters can be removed
    elimination_done = false;
    elim_iter = 0;
    
    fprintf('    > Starting Backward Elimination...\n');
    
    while ~elimination_done
        elim_iter = elim_iter + 1;
        
        % 3.3.1) Select columns for currently active parameters
        A_curr = A_full(:, mask_params);
        
        % 3.3.2) Solve LS: (A' * Q^-1 * A) * x = A' * Q^-1 * y
        W_A = Q_ll \ A_curr; 
        N   = A_curr' * W_A;
        rhs = W_A' * y0_lsc;
        
        % Estimate x (with stability check)
        if cond(N) > 1e12
             N_inv = pinv(N);
             x_curr = N_inv * rhs;
        else
             N_inv = inv(N);
             x_curr = N \ rhs;
        end
        
        % 3.3.3) Compute statistics
        % Residuals
        v = y0_lsc - A_curr * x_curr;
        
        % Weighted Sum of Squared Residuals (WSSR) (v' * inv(Q) * v)
        W_v = Q_ll \ v;
        WSSR = v' * W_v;
        
        % Degrees of Freedom
        n_params_curr = sum(mask_params);
        dof = length(y0_lsc) - n_params_curr;
        
        % A Posteriori Variance Factor
        sigma0_sq = WSSR / dof;
        
        % Covariance of Parameters (Cxx)
        C_xx = sigma0_sq * N_inv;
        
        % T-Statistics
        x_std = sqrt(diag(C_xx));
        t_stats = x_curr ./ x_std;
        
        % Critical Value (two-tailed)
        try
            t_crit = tinv(1 - alpha_sig/2, dof);
        catch
            t_crit = 1.96; 
        end
        
        % 3.3.4) Identify the smallest one
        [min_t, min_idx_local] = min(abs(t_stats));
        
        % Map local index back to global index to identify what it is
        active_global_indices = find(mask_params);
        worst_param_global_idx = active_global_indices(min_idx_local);
        
        % 3.3.5) Decision Logic
        if min_t < t_crit
            % Determine what is removed for the log
            if worst_param_global_idx <= n_trend_cols
                param_type = 'TREND';
                param_desc = sprintf('Coeff %d', worst_param_global_idx);
            else
                param_type = 'MEAN';
                % Map to epoch index
                epoch_idx = worst_param_global_idx - n_trend_cols;
                param_desc = sprintf('Epoch %d', epoch_idx);
            end
            
            fprintf('      [Elim %d] Removing %s (%s). t=%.3f < %.3f\n', ...
                elim_iter, param_type, param_desc, min_t, t_crit);
            
            % Remove the parameter
            mask_params(worst_param_global_idx) = false;
            
            % Safety check: if everything removed, stop
            if ~any(mask_params)
                elimination_done = true;
            end
        else
            % All remaining parameters are significant
            elimination_done = true;
        end
    end
        
    
    % 3.4) Reconstruct full vectors for the rest of the processing
    x_est_all = zeros(size(A_full, 2), 1);
    x_est_all(mask_params) = x_curr;
    
    switch detrend_method
        case 'residualAtm'
        % 3.4.1) Extract Global Orbit Trend
        x_trend = x_est_all(1:n_trend_cols);
        
        % Reconstruct Global Trend Grid
        A_red_trend = s_full_centered;
        trend_full_col = A_red_trend * x_trend;
        trend_grid_full = reshape(trend_full_col, length(t_relIN), [])';
    
        % 3.4.2) Extract Atmosphere
        x_atmos_all = x_est_all(n_trend_cols+1:end);
        
        % Reconstruct Atmosphere Grid
        if ~use_inclined_means
            % Case: Constant means
            % x_atmos_all is [n_epochs x 1]
            spatialMeans_grid = repmat(x_atmos_all', size(residuals_grid_full, 1), 1);
            
        else
            % Case: Inclined Planes (Ramps)
            % x_atmos_all is [2*n_epochs x 1]
            % Reshape to [2 x n_epochs] -> Row 1: Bias, Row 2: SlopeS
            atmos_params = reshape(x_atmos_all, 2, n_epochs);
            
            biases  = atmos_params(1, :); % [1 x n_epochs]
            slopesS = atmos_params(2, :);
            
            % Construct the full grid correction
            % Atmos(t,x,y) = Bias(t) + Ss(t)*s
            
            % Grid S (centered)
            S_grid = reshape(s_full_centered, length(t_relIN), [])'; % [n_PS x n_epochs]
            
            % Expand parameters to grid size
            B_grid = repmat(biases, size(residuals_grid_full, 1), 1);
            Ss_grid = repmat(slopesS, size(residuals_grid_full, 1), 1);
            
            spatialMeans_grid = B_grid + (Ss_grid .* S_grid);
        end
    
    end

end


% 3.5) Polynomial surface fitting (empirical)
% If 'residualAtm' was run, the polynomial is fitted on the residuals of that model.
% If 'cleanObs' was selected, it is fitted on the raw y0_lsc.
if exist('x_est_all', 'var') && exist('A_full', 'var')
    % remove the previously estimated Orbit + Atmosphere
    y_for_poly = y0_lsc - A_full * x_est_all;
else
    % no previous subtraction
    y_for_poly = y0_lsc;
end

% 3.5.1) Build design matrix A_poly
% (using the vectors corresponding to current VALID observations)
t_n = t_centered;
s_n = s_centered;

% Start with intercept
A_poly = ones(length(y0_lsc), 1);
param_names = {'Bias'};

% Degree 1
if poly_degree >= 1
    A_poly = [A_poly, s_n, t_n];
    param_names = [param_names, {'s', 't'}];
end

% Degree 2
if poly_degree >= 2
    A_poly = [A_poly, s_n.^2, t_n.^2, s_n.*t_n];
    param_names = [param_names, {'s^2', 't^2', 'st'}];
end

% Degree 3
if poly_degree >= 3
    A_poly = [A_poly, s_n.^3, t_n.^3, ...
              s_n.^2.*t_n, t_n.^2.*s_n];
    param_names = [param_names, {'s^3', 't^3', 's^2t', 'st^2'}];
end

% 3.5.2) Backward elimination for polynomial
mask_poly = true(size(A_poly, 2), 1);
elimination_done = false;
elim_iter = 0;

fprintf('    > Starting Backward Elimination for Polynomial...\n');
        
while ~elimination_done
    elim_iter = elim_iter + 1;
    
    % Select active columns
    A_curr = A_poly(:, mask_poly);
    
    % Solve LS
    W_A = Q_ll \ A_curr; 
    N   = A_curr' * W_A;
    rhs = W_A' * y_for_poly;
    
    if cond(N) > 1e12
         N_inv = pinv(N);
         x_curr = N_inv * rhs;
    else
         N_inv = inv(N);
         x_curr = N \ rhs;
    end
    
    % Statistics
    v = y_for_poly - A_curr * x_curr;
    WSSR = v' * (Q_ll \ v);
    dof = length(y_for_poly) - sum(mask_poly);
    sigma0_sq = WSSR / dof;
    C_xx = sigma0_sq * N_inv;
    t_stats = x_curr ./ sqrt(diag(C_xx));
    
    try
        t_crit = tinv(1 - alpha_sig/2, dof);
    catch
        t_crit = 1.96; 
    end
    
    % Identify worst parameter
    [min_t, min_idx_local] = min(abs(t_stats));
    active_indices = find(mask_poly);
    worst_global_idx = active_indices(min_idx_local);
    
    % Decision
    if min_t < t_crit
        fprintf('      [Poly Elim %d] Removing %s. t=%.3f < %.3f\n', ...
            elim_iter, param_names{worst_global_idx}, min_t, t_crit);
        mask_poly(worst_global_idx) = false;
        
        if ~any(mask_poly), elimination_done = true; end
    else
        elimination_done = true;
    end
end

% 3.5.3) Store coefficients and reconstruct full grid
% Store the final coefficients (zeros for removed params)
x_poly_coeffs = zeros(size(A_poly, 2), 1);
x_poly_coeffs(mask_poly) = x_curr;

% Reconstruct on the full grid (using s_full_centered, t_full_centered)
np_full = size(residuals_grid_full, 1);
nt_full = size(residuals_grid_full, 2);

% Flatten full grid coords (ensuring alignment with how A_poly was built)
s_f = s_full_centered;   % [np*nt x 1]
t_f = t_full_centered;   % [np*nt x 1]

% Build A_poly_full
A_poly_full = ones(length(s_f), 1);

if poly_degree >= 1
    A_poly_full = [A_poly_full, s_f, t_f];
end

if poly_degree >= 2
    A_poly_full = [A_poly_full, s_f.^2, t_f.^2, s_f.*t_f];
end

if poly_degree >= 3
    A_poly_full = [A_poly_full, s_f.^3, t_f.^3, s_f.^2.*t_f, t_f.^2.*s_f];
end

% Calculate trend on full grid
trend_poly_col = A_poly_full * x_poly_coeffs;
trend_poly_grid = reshape(trend_poly_col, nt_full, np_full)';   % reshape back to [np x nt]

% 3.6) Final residuals (signal)
% The "clean" signal is:
% original - (orbit+atmos if exists) - (polynomial)

if exist('trend_grid_full', 'var') && exist('spatialMeans_grid', 'var')
    % If atmos method was used:
    % total model = orbit + atmos + poly
    full_model_grid = trend_grid_full + spatialMeans_grid + trend_poly_grid;
else
    % If only poly was used:
    full_model_grid = trend_poly_grid;
end

% Final cleaned time series for splies
data_final_clean = residuals_grid_full - full_model_grid;


disp('Least Squares model reduction completed.');




%% 4) SPLINES INTERPOLATION
disp('======== Step 4 ========');
disp('Splines interpolation started...');

% 4.1) Define noise variance
switch varNoise_method
    case 'auto'
        % METHOD: Pairwise Difference (Allan Variance approximation)
        % It assumes the unmodeled signal is temporally correlated (smooth),
        % while the noise is white (uncorrelated).
        % Taking the difference between consecutive residuals cancels the signal
        % and leaves 2*noise_variance.
        
        % Compute differences along the time dimension
        diffs = diff(res_spl, 1, 2); 
        
        % Calculate variance of these differences
        var_diffs = var(diffs(:), 'omitnan');
        
        % The noise variance is half of the difference variance
        varNoise = var_diffs / 2;
    case 'manual'
        varNoise = varNoise_manual;
end


% 4.2) Define domain ranges
x_coords_sort = grid_s_flt_col;
y_coords_sort = grid_t_flt_col;
data_final_clean_col = reshape(data_final_clean', [], 1);
data_final_clean_col(isnan(data_final_clean_col)) = [];

x_min = min(x_coords_sort);
x_max = max(x_coords_sort);
y_min = min(y_coords_sort);
y_max = max(y_coords_sort);
x_range = x_max - x_min;
y_range = y_max - y_min;

% Compute average distance between PS
avg_dist_s = median(diff(s_SAR));
avg_dist_t = median(diff(t_relIN));

% 20% of domain or 20x median spacing
range_s = max(x_range / 5, avg_dist_s * 20);  
range_t = max(y_range / 5, avg_dist_t * 20);

% Enforce minimum ranges
range_s = max(range_s, min(x_range / 10, 2*round(mode(diff(s_SAR)))));
range_t = max(range_t, min(y_range / 10, 2*round(mode(diff(t_relIN)))));


% 4.3) Processing domains
% Initial estimate for num_row and num_col based on correlation ranges
base_num_row = max(2, round(x_range / range_s));  % knots spaced at ~range_s
base_num_col = max(2, round(y_range / range_t));  % knots spaced at ~range_t

% Define min and max number of splines to test in each direction
min_n_row = max(2, round(base_num_row * 0.5));  % 50% of base
max_n_row = min(round(base_num_row * 2), round(length(s_SAR) / 4));  % 200% or 1/4 points
min_n_col = max(2, round(base_num_col * 0.5));
max_n_col = min(round(base_num_col * 2), round(length(t_relIN) / 4));

% Compute dynamic step sizes based on range width
row_range = max_n_row - min_n_row;
col_range = max_n_col - min_n_col;
row_step = max(1, ceil(row_range / 6));
col_step = max(1, ceil(col_range / 6));

% Ensure base_num_row and base_num_col are included
row_test = unique([min_n_row:row_step:max_n_row, base_num_row]);
col_test = unique([min_n_col:col_step:max_n_col, base_num_col]);

% Cap at 7 values per dimension to limit combinations
if length(row_test) > 7
    row_test = row_test(round(linspace(1, length(row_test), 7)));
end
if length(col_test) > 7
    col_test = col_test(round(linspace(1, length(col_test), 7)));
end


% 4.4) Write the input file
gS_filename = 'DET1D';
writematrix([x_coords_sort, y_coords_sort, data_final_clean_col], strcat(gS_input_path, filesep, gS_filename, '.txt'), 'Delimiter', 'space')

% Define processing variables
ext = 'bic';                                    % ext variable based on splines type

% Job file creation
data_dim = 2;                                   % Dimension of dataset (1D)
type_spl = 2;                                   % Type of splines (bicubic)
file_inp = strcat(gS_filename, '.txt');         % Input filename (with extension)
num_obs  = length(data_final_clean_col);                    % Number of observations (PS)
x_min;                                          % First x-coordinate
x_max;                                          % Last x-coordinate
y_min;                                          % First y-coordinate
y_max;                                          % Last y-coordinate
num_sig  = 10;                                  % Number of significant digits
    
% Validate manual inputs
switch num_spl_method
    case 'manual'
    if ~isnumeric(num_spl_row_manual) || ~isscalar(num_spl_row_manual) || num_spl_row_manual < min_n_row || mod(num_spl_row_manual, 1) ~= 0
        warning('Invalid num_spl_row_manual (%s). Reverting to auto.', num2str(num_spl_row_manual));
        num_spl_method = 'auto';
    end
    if ~isnumeric(num_spl_col_manual) || ~isscalar(num_spl_col_manual) || num_spl_col_manual < min_n_col || mod(num_spl_col_manual, 1) ~= 0
        warning('Invalid num_spl_col_manual (%s). Reverting to auto.', num2str(num_spl_col_manual));
        num_spl_method = 'auto';
    end
end
switch lambda_method
    case 'manual'
    if ~isnumeric(lambda_manual) || ~isscalar(lambda_manual) || lambda_manual <= 0
        warning('Invalid lambda_manual (%s). Reverting to auto.', num2str(lambda_manual));
        lambda_method = 'auto';
    end
end


% 4.5) Initialize variables for candidate solution
switch num_spl_method
    case 'manual'
        num_row_candidates = 1;
        num_col_candidates = 1;
    case 'auto'
        num_row_candidates = length(row_test);
        num_col_candidates = length(col_test);
end
num_candidates = num_row_candidates * num_col_candidates;
PS_data_spline_tmp = NaN(num_obs, num_candidates);
s02_spl_tmp = NaN(1, num_candidates);
var_resSpl = NaN(num_candidates, 1);
MDL_tmp = NaN(num_candidates, 1);
F_pval_tmp = NaN(num_candidates, 1);
chi2_pval_tmp = NaN(num_candidates, 1);
lambda_it = NaN(num_candidates, 1);
std_spl_gS = cell(num_candidates, 1);
row_col_pairs = NaN(num_candidates, 2);
knots_x_tmp = cell(num_candidates, 1);
knots_y_tmp = cell(num_candidates, 1);
raster_tmp = cell(num_candidates, 1);

% Determine spline combinations to test
switch num_spl_method
    case 'manual'
        rows_to_test = num_spl_row_manual;
        cols_to_test = num_spl_col_manual;
    case 'auto'
        rows_to_test = row_test;
        cols_to_test = col_test;
end


% 4.6) Iterate over num_row and num_col combinations
kk = 1;
for r = rows_to_test
    for c = cols_to_test

        % Incremental definition of spline number
        num_row = r;
        num_col = c;

        % Number of parameters for bicubic splines
        n_params = (num_row + 2) * (num_col + 2);

        % Check if num_obs is sufficient
        if num_obs <= n_params
            switch num_spl_method
                case 'auto'
                    warning('num_obs (%d) <= n_params (%d), iteration %d. Skipping.', ...
                        num_obs, n_params, kk);
                case 'manual'
                    error('num_obs (%d) <= n_params (%d). Reduce the number of splines.', ...
                        num_obs, n_params);
            end
            continue;
        end

        % Define output filename for this combination
        file_out = sprintf('%s_%s_r%d_c%d', gS_filename, ext, num_row, num_col);

        % Regularization parameter (λ)
        switch lambda_method
            case 'auto'
            tauGrid_x = x_range / (num_row - 1);
            tauGrid_y = y_range / (num_col - 1);
            deltaGrid = mean([tauGrid_x, tauGrid_y]);
            lambda = lambdaSplines2D([x_coords_sort, y_coords_sort, data_final_clean_col], ...
                                     varNoise, deltaGrid, type_spl);
            case 'manual'
            lambda = lambda_manual;
        end

        % Write the job file
        jobFile_analysis2D(data_dim, type_spl, file_inp, file_out, num_obs, ...
                           num_row, num_col, x_min, x_max, y_min, y_max, ...
                           lambda, num_sig, gS_input_path, gS_output_path, gS_job_path);

        % Run geoSplinter_analysis
        job_execution = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_analysis'), ...
                                fullfile('.', gS_job_path, [file_out, '.job']));
        status = system(job_execution);
        if status ~= 0
            error('Error executing geoSplinter_analysis for file: %s', file_out);
        end

        % Import results
            %[~, ~, ~,~] = geoSplinter(fullfile('.', gS_output_path, file_out), ext);
        [gS_out_data_tmp, gS_out_raster, ~, ~, gS_out_nodes] = geoSplinter_noFig(...
            fullfile('.', gS_output_path, file_out), ext);

        % Import knots
        knots_x = gS_out_nodes(1,:,1);
        knots_y = gS_out_nodes(:,1,2)';
        deltaX = knots_x(2) - knots_x(1);
        deltaY = knots_y(1) - knots_y(2);

        % Define support functions for geoSplinter
        phi_3 = @(csi) ((2 - csi).^3 - 4 * (1 - csi).^3) / 6;
        phi_4 = @(csi) (2 - csi).^3 / 6;
        phi_33 = @(csi_x, csi_y) phi_3(csi_x) .* phi_3(csi_y);
        phi_34 = @(csi_x, csi_y) phi_3(csi_x) .* phi_4(csi_y);
        phi_43 = @(csi_x, csi_y) phi_4(csi_x) .* phi_3(csi_y);
        phi_44 = @(csi_x, csi_y) phi_4(csi_x) .* phi_4(csi_y);

        % Initialize full design matrix
        A_spline = zeros(length(x_coords_sort), n_params);
        
        % Compute design matrix
        for n = 1:length(x_coords_sort)
            x = x_coords_sort(n);
            y = y_coords_sort(n);
            
            % Find knot interval for x: knots_x(i_x) <= x < knots_x(i_x+1)
            i_x = find(x >= knots_x(1:end-1) & x < knots_x(2:end), 1);
            if isempty(i_x)
                if x <= knots_x(1)
                    i_x = 1;
                elseif x >= knots_x(end)
                    i_x = length(knots_x) - 1;
                end
            else
                i_x = i_x - 1; 
            end
            
            % Find knot interval for y: knots_y(i_y) >= y > knots_y(i_y+1) (decreasing)
            i_y = find(y <= knots_y(1:end-1) & y > knots_y(2:end), 1);
            if isempty(i_y)
                if y >= knots_y(1)
                    i_y = 1;
                elseif y <= knots_y(end)
                    i_y = length(knots_y) - 1;
                end
            else
                i_y = i_y - 1; 
            end
            
            % Check valid indices (mimicking C++: i_x >= -2 && i_x <= xNum && i_y >= -2 && i_y <= yNum)
            if (i_x >= -2 && i_x <= length(knots_x) && i_y >= -2 && i_y <= length(knots_y))
                % Compute local coordinates
                csi_x = (x - knots_x(i_x + 1)) / deltaX;
                csi_y = (y - knots_y(i_y + 1)) / deltaY;
                
                % Compute alpha coefficients (4x4)
                alpha = zeros(4, 4);
                alpha(1,1) = phi_44(1 + csi_x, 1 + csi_y);
                alpha(1,2) = phi_43(1 + csi_x, csi_y);
                alpha(1,3) = phi_43(1 + csi_x, 1 - csi_y);
                alpha(1,4) = phi_44(1 + csi_x, 2 - csi_y);
                alpha(2,1) = phi_34(csi_x, 1 + csi_y);
                alpha(2,2) = phi_33(csi_x, csi_y);
                alpha(2,3) = phi_33(csi_x, 1 - csi_y);
                alpha(2,4) = phi_34(csi_x, 2 - csi_y);
                alpha(3,1) = phi_34(1 - csi_x, 1 + csi_y);
                alpha(3,2) = phi_33(1 - csi_x, csi_y);
                alpha(3,3) = phi_33(1 - csi_x, 1 - csi_y);
                alpha(3,4) = phi_34(1 - csi_x, 2 - csi_y);
                alpha(4,1) = phi_44(2 - csi_x, 1 + csi_y);
                alpha(4,2) = phi_43(2 - csi_x, csi_y);
                alpha(4,3) = phi_43(2 - csi_x, 1 - csi_y);
                alpha(4,4) = phi_44(2 - csi_x, 2 - csi_y);
                
                % Assign to design matrix
                for k = -1:2
                    for h = -1:2
                        if (i_x + k >= 0 && i_x + k <= length(knots_x) - 1 && ...
                            i_y + h >= 0 && i_y + h <= length(knots_y) - 1)
                            % Column index (column-major, matching order(i_x, i_y) = i_y + i_x * yNum)
                            col = (i_y + h) + (i_x + k) * length(knots_y) + 1;
                            A_spline(n, col) = alpha(k + 2, h + 2);
                        end
                    end
                end
            end
        end

        % Compute the normal matrix
        N_spline = A_spline' * A_spline;

        % Compute residuals and variance
        estimated_obs = gS_out_data_tmp(:, 4);
        residuals = data_final_clean_col - estimated_obs;
        SSR = sum(residuals.^2);
        s02_est = SSR / (num_obs - n_params);
        Cyy_est = s02_est * (A_spline * (N_spline \ A_spline'));

        % Store results
        PS_data_spline_tmp(:, kk) = estimated_obs;
        s02_spl_tmp(kk) = s02_est;
        var_resSpl(kk) = s02_est;
        lambda_it(kk) = lambda;
        row_col_pairs(kk, :) = [num_row, num_col];
        knots_x_tmp{kk} = knots_x;
        knots_y_tmp{kk} = knots_y;
        raster_tmp{kk} = gS_out_raster;
        std_spl_gS{kk} = sqrt(diag(Cyy_est));

        % Compute selection metrics for 'auto' method
        if strcmp(num_spl_method, 'auto')
            % MDL index
            MDL_tmp(kk) = log(num_obs) * (num_row * num_col) / num_obs + log(max(SSR, eps));

            % Chi-squared test p-value
            if num_obs > n_params && varNoise > 0
                chi2_stat = SSR / varNoise;
                dof = num_obs - n_params;
                chi2_pval_tmp(kk) = 1 - chi2cdf(chi2_stat, dof);
            else
                chi2_pval_tmp(kk) = NaN;
            end
        end

        kk = kk + 1;
    end
end

% Compute F-test p-values
% (compare num_row/num_col vs. num_row_prev/num_col_prev)
switch num_spl_method
    case 'auto'
        % Identify the unique grid steps actually present
        u_rows = unique(row_col_pairs(:, 1));
        u_cols = unique(row_col_pairs(:, 2));

        for kk = 1:num_candidates
            % Current model parameters
            curr_r = row_col_pairs(kk, 1);
            curr_c = row_col_pairs(kk, 2);
            
            % Need to find the "simpler model" (predecessor)
            % Find where current row/col sit in the unique list
            r_idx_list = find(u_rows == curr_r);
            c_idx_list = find(u_cols == curr_c);
            
            % Can only compare if there is a previous step in both dimensions
            if r_idx_list > 1 && c_idx_list > 1
                % Get the actual previous values
                prev_r = u_rows(r_idx_list - 1);
                prev_c = u_cols(c_idx_list - 1);
                
                % Find the index of this simpler model in the main list
                simpler_idx = find(row_col_pairs(:,1) == prev_r & row_col_pairs(:,2) == prev_c, 1);
                
                % Check if simpler_idx exists and if it was successfully computed (not skipped)
                if ~isempty(simpler_idx) && ~isnan(s02_spl_tmp(simpler_idx)) && ~isnan(s02_spl_tmp(kk))
                    
                    % Complex model (current)
                    n_params_complex = (curr_r + 2) * (curr_c + 2);
                    dof_complex = num_obs - n_params_complex;
                    SSR_complex = s02_spl_tmp(kk) * dof_complex;
                    
                    % Simple model (predecessor)
                    n_params_simple = (prev_r + 2) * (prev_c + 2);
                    dof_simple = num_obs - n_params_simple;
                    SSR_simple = s02_spl_tmp(simpler_idx) * dof_simple;
                    
                    % F-Statistic
                    % (Reduction in error / Cost of new params) / Error of complex
                    numerator = (SSR_simple - SSR_complex) / (dof_simple - dof_complex);
                    denominator = s02_spl_tmp(kk);   % MSE of complex
                    
                    if dof_simple > dof_complex && denominator > 0
                        F_stat = numerator / denominator;
                        % P-value: probability that improvement is just random noise
                        F_pval_tmp(kk) = 1 - fcdf(F_stat, dof_simple - dof_complex, dof_complex);
                    else
                        F_pval_tmp(kk) = NaN;
                    end
                else
                    F_pval_tmp(kk) = NaN;   % cannot compare (one model missing/failed)
                end
            else
                % This is the simplest model (base case), cannot test downwards
                F_pval_tmp(kk) = NaN; 
            end
        end
    
        % Set p-value for the absolute simplest model to 0 (so we can pick it if all else fails)
        base_idx = find(row_col_pairs(:,1) == u_rows(1) & row_col_pairs(:,2) == u_cols(1));
        if ~isempty(base_idx)
            F_pval_tmp(base_idx) = 0; 
        end
end

% 4.7) Select number of splines
switch num_spl_method
    case 'manual'
        num_row = num_spl_row_manual;
        num_col = num_spl_col_manual;
        var_idx = find(row_col_pairs(:,1) == num_row & row_col_pairs(:,2) == num_col, 1);
        if isempty(var_idx)
            error('Manual num_row=%d, num_col=%d not found in tested combinations.', num_row, num_col);
        end
    case 'auto'
        % Identify which models were actually computed (and not skipped/failed)
        valid_computed_mask = ~cellfun(@isempty, std_spl_gS);
        valid_indices = find(valid_computed_mask);
        if isempty(valid_indices)
            error('ALL spline models were skipped! Your num_obs (%d) is likely too small for the requested spline density.', num_obs);
        end
        switch spline_method
            case 'variance'
                % Heuristic: pick model with variance close to 80th percentile (avoid overfitting)
                % Only look at valid computed variances
                valid_vars = var_resSpl(valid_indices);
                var_thrs = prctile(valid_vars, 80);
                diffs = valid_vars - var_thrs;
                [~, min_diff_idx] = min(diffs);
                if isempty(min_diff_idx) || isinf(diffs(min_diff_idx))
                    var_idx = valid_indices(1); 
                else
                    var_idx = valid_indices(min_diff_idx);
                end
            case 'MDL'
                % MDL: find the minimum MDL among computed models
                % Initialize with Inf so skipped models (NaN) are ignored
                MDL_clean = inf(size(MDL_tmp));
                MDL_clean(valid_indices) = MDL_tmp(valid_indices);
                [~, var_idx] = min(MDL_clean);
            case 'F_test'
                % F-test: find the most complex model that is statistically significant
                valid_pvals = F_pval_tmp;
                valid_pvals(~valid_computed_mask) = NaN;
                if all(isnan(valid_pvals))
                    warning('F-test: No valid comparisons found. Defaulting to simplest computed model.');
                    var_idx = find(valid_computed_mask, 1, 'first');
                else
                    candidate_idxs = find(valid_pvals < alpha_sig);
                    if isempty(candidate_idxs)
                         var_idx = find(valid_computed_mask, 1, 'first');
                    else
                         var_idx = max(candidate_idxs);
                    end
                end
            case 'chi2_test'
                % Chi2_test: check if p-values fall within the acceptance range [alpha, 1-alpha]
                valid_chi2 = chi2_pval_tmp;
                valid_chi2(~valid_computed_mask) = NaN;
                idx = find(valid_chi2 >= alpha_sig/2 & valid_chi2 <= (1 - alpha_sig/2), 1, 'first');
                if isempty(idx)
                    warning('Chi2 test failed to find a match. Defaulting to simplest valid model.');
                    var_idx = valid_indices(1); 
                else
                    var_idx = idx;
                end
            otherwise
                error('Invalid spline_method: %s. Choose ''variance'', ''MDL'', ''F_test'', ''chi2_test''.', spline_method);
        end
        num_row = row_col_pairs(var_idx, 1);
        num_col = row_col_pairs(var_idx, 2);
end

% Get the result from var_idx
PS_data_spline = PS_data_spline_tmp(:, var_idx);
res_spl = data_final_clean_col - PS_data_spline;
lambda_spl = lambda_it(var_idx);
std_spl = std_spl_gS{var_idx}; 
knots_x = knots_x_tmp{var_idx};
knots_y = knots_y_tmp{var_idx};
raster = raster_tmp{var_idx};
s02_spl = var_resSpl(var_idx);


% 4.8) Splines solution estimation on the full query grid
% create the full grid of estimation distances and times
t_full = (t_relIN(1):step_t:t_relIN(end))'; % [n_t x 1]
s_full = s_full_SAR; 

% Define dimensions
n_t_query = length(t_full);
n_s_query = length(s_full);

% Use ndgrid to create Time-Major matrices
% T_grid will vary down the rows (Time-Major)
[T_grid, S_grid] = ndgrid(t_full, s_full);

% Flatten to column vectors
% The (:) operator flattens column-by-column. 
% Since ndgrid puts T on rows, T varies fastest.
grid_t_full_col = T_grid(:);
grid_s_full_col = S_grid(:);

clear T_grid S_grid;

% Compute the solution on the full grid
xy_est = [grid_s_full_col, grid_t_full_col];
file_out = sprintf('%s_bic_r%d_c%d', gS_filename, num_row, num_col);

% Synthesize splines solution on xy_est
gS_est = 'st_est';
writematrix(xy_est, fullfile(gS_input_path, [gS_est, '.txt']), 'Delimiter', 'space');

% Job file creation for synthesis
file_spl = file_out;
file_est = strcat(gS_est, '.txt');
file_syn = sprintf('%s_bic_est', gS_filename);
jobFile_synthesis(data_dim, type_spl, file_spl, file_est, file_syn, gS_input_path, gS_output_path, gS_job_path, gS_synth_path);

% Run geoSplinter_synthesis
job_execution_syn = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_synthesis'), fullfile('.', gS_job_path, [file_syn, '.job']));
status = system(job_execution_syn);
if status ~= 0
    error('Error executing geoSplinter_synthesis for file: %s', file_syn);
end

% Results import
spl_est_full = readmatrix(fullfile(gS_synth_path, file_syn));
spl_est = spl_est_full(:, 3);

% Reshape to a grid
spl_est_grid = reshape(spl_est, length(t_full), [])';
figure; imagesc(spl_est_grid)

% 4.9) Restore the removed polynomial surface
t_est = xy_est(:, 2);
s_est_vec = xy_est(:, 1);

% Centering (must match the loop logic)
s_est_centered = s_est_vec - s0;
t_est_centered = t_est - t0;

% Estimate the polynomial at query grid
t_n = t_est_centered;
s_n = s_est_centered;

A_poly = ones(length(t_n), 1);

if poly_degree >= 1
    A_poly = [A_poly, s_n, t_n];
end

if poly_degree >= 2
    A_poly = [A_poly, s_n.^2, t_n.^2, s_n.*t_n];
end

if poly_degree >= 3
    A_poly = [A_poly, s_n.^3, t_n.^3, ...
              s_n.^2.*t_n, t_n.^2.*s_n];
end

z_poly_est = A_poly * x_poly_coeffs;

% Reshape trend to grid
poly_fit_full = reshape(z_poly_est, n_t_query, n_s_query)';


% 4.10) Get the final signal
final_signal = spl_est_grid + poly_fit_full;


% 4.11) Apply reference shift (relative to first epoch, t=0)
% This ensures the deformation starts at 0 for all pixels.
final_signal_shift = final_signal - final_signal(:,1);

% Also shift the components for plotting consistency
spl_est_grid_shift = spl_est_grid - spl_est_grid(:,1);
poly_fit_full_grid_shift = poly_fit_full - poly_fit_full(:,1);


% 4.12) Prediction on observation points
% Define target coordinates (full original grid)
% These variables were created in Step 2
target_t = t_full_0; % [np*nt x 1]
target_s = s_full_0;

% 4.12.1) Compute splines component on full grid
xy_obs_target = [target_s, target_t];

% Write target coordinates to file
gS_est_obs = 'st_obs_full';
writematrix(xy_obs_target, fullfile(gS_input_path, [gS_est_obs, '.txt']), 'Delimiter', 'space');

% Job file creation for synthesis
file_syn_obs = sprintf('%s_bic_est_obs', gS_filename);
jobFile_synthesis(data_dim, type_spl, file_spl, [gS_est_obs, '.txt'], file_syn_obs, ...
                  gS_input_path, gS_output_path, gS_job_path, gS_synth_path);

% Run geoSplinter_synthesis
job_execution_syn_obs = sprintf('%s < %s', fullfile(gS_dir, 'geoSplinter_synthesis'), ...
                                fullfile('.', gS_job_path, [file_syn_obs, '.job']));
status = system(job_execution_syn_obs);
if status ~= 0
    error('Error executing geoSplinter_synthesis for Obs Grid: %s', file_syn_obs);
end

% Import results
spl_est_obs_full = readmatrix(fullfile(gS_synth_path, file_syn_obs));
signal_spl_full = spl_est_obs_full(:, 3);

% 4.12.2) Compute polynomial component at observation points
% centering (as Step 3)
t_target_n = target_t - t0;
s_target_n = target_s - s0;

% Build design matrix (A_poly_target)
A_poly_target = ones(length(t_target_n), 1);

if poly_degree >= 1
    A_poly_target = [A_poly_target, s_target_n, t_target_n];
end

if poly_degree >= 2
    A_poly_target = [A_poly_target, s_target_n.^2, t_target_n.^2, ...
                  s_target_n.*t_target_n];
end
if poly_degree >= 3
    A_poly_target = [A_poly_target, s_target_n.^3, t_target_n.^3, ...
                  s_target_n.^2.*t_target_n, t_target_n.^2.*s_target_n];
end

% Apply coefficients (x_poly_coeffs from Step 3)
signal_poly_full = A_poly_target * x_poly_coeffs;

% 4.12.3) Final summation and reshape
final_signal_orig_col = signal_spl_full + signal_poly_full;
final_signal_orig = reshape(final_signal_orig_col, nt, np)';

disp('Splines interpolation completed.');




%% 5) TOTAL VARIANCE PROPAGATION
disp('======== Step 5 ========');
disp('Variance propagation started...');

% 5.1) Setup target coordinates
% xy_est is [n_total x 2] -> [Time, S]
n_2D_total = size(xy_est, 1); 


% 5.2) Construct design matrices (A_est_2D and A_obs)
% Coordinates for observations (centered)
s_obs_n = grid_s_flt_col - s0;
t_obs_n = grid_t_flt_col - t0;

% Build Full A_poly_obs
A_poly_obs = ones(length(s_obs_n), 1);
if poly_degree >= 1
    A_poly_obs = [A_poly_obs, s_obs_n, t_obs_n];
end
if poly_degree >= 2
    A_poly_obs = [A_poly_obs, s_obs_n.^2, t_obs_n.^2, s_obs_n.*t_obs_n];
end
if poly_degree >= 3
    A_poly_obs = [A_poly_obs, s_obs_n.^3, t_obs_n.^3, ...
                  s_obs_n.^2.*t_obs_n, t_obs_n.^2.*s_obs_n];
end

% Apply the mask from Backward Elimination
A_poly_obs_red = A_poly_obs(:, mask_poly);

% Compute Polynomial Parameter Covariance Matrix (C_xx_poly)
% N = A' * W * A
% C_xx = sigma0^2 * inv(N)
W_A_poly = Q_ll \ A_poly_obs_red;
N_poly   = A_poly_obs_red' * W_A_poly;
if cond(N_poly) > 1e12
    C_xx_poly = sigma0_sq * pinv(N_poly);
else
    C_xx_poly = sigma0_sq * inv(N_poly);
end

% Reconstruct A_poly_est (Design matrix on Prediction Grid)
% (already built in Step 4.9 as A_poly, but need to apply the mask)
A_poly_est_red = A_poly(:, mask_poly);

% Propagate to grid: diag(A * C_xx * A')
var_poly_grid = sum((A_poly_est_red * C_xx_poly) .* A_poly_est_red, 2);


% 5.3) Variance of the spline component
% 5.3.1) Reconstruct N_spline (observation normal matrix)
% Rely on the knots and sigma0 (s02_spl) from the best model in Step 4

% Parameters from best fit
n_params_spl = (num_row + 2) * (num_col + 2);
deltaX = knots_x(2) - knots_x(1);
deltaY = knots_y(1) - knots_y(2);

    A_spline_obs = zeros(length(x_coords_sort), n_params_spl);
        
    % Compute design matrix
    for n = 1:length(x_coords_sort)
        x = x_coords_sort(n);
        y = y_coords_sort(n);
        
        % Find knot interval for x: knots_x(i_x) <= x < knots_x(i_x+1)
        i_x = find(x >= knots_x(1:end-1) & x < knots_x(2:end), 1);
        if isempty(i_x)
            if x <= knots_x(1)
                i_x = 1;
            elseif x >= knots_x(end)
                i_x = length(knots_x) - 1;
            end
        else
            i_x = i_x - 1; 
        end
        
        % Find knot interval for y: knots_y(i_y) >= y > knots_y(i_y+1) (decreasing)
        i_y = find(y <= knots_y(1:end-1) & y > knots_y(2:end), 1);
        if isempty(i_y)
            if y >= knots_y(1)
                i_y = 1;
            elseif y <= knots_y(end)
                i_y = length(knots_y) - 1;
            end
        else
            i_y = i_y - 1; 
        end
        
        % Check valid indices (mimicking C++: i_x >= -2 && i_x <= xNum && i_y >= -2 && i_y <= yNum)
        if (i_x >= -2 && i_x <= length(knots_x) && i_y >= -2 && i_y <= length(knots_y))
            % Compute local coordinates
            csi_x = (x - knots_x(i_x + 1)) / deltaX;
            csi_y = (y - knots_y(i_y + 1)) / deltaY;
            
            % Compute alpha coefficients (4x4)
            alpha = zeros(4, 4);
            alpha(1,1) = phi_44(1 + csi_x, 1 + csi_y);
            alpha(1,2) = phi_43(1 + csi_x, csi_y);
            alpha(1,3) = phi_43(1 + csi_x, 1 - csi_y);
            alpha(1,4) = phi_44(1 + csi_x, 2 - csi_y);
            alpha(2,1) = phi_34(csi_x, 1 + csi_y);
            alpha(2,2) = phi_33(csi_x, csi_y);
            alpha(2,3) = phi_33(csi_x, 1 - csi_y);
            alpha(2,4) = phi_34(csi_x, 2 - csi_y);
            alpha(3,1) = phi_34(1 - csi_x, 1 + csi_y);
            alpha(3,2) = phi_33(1 - csi_x, csi_y);
            alpha(3,3) = phi_33(1 - csi_x, 1 - csi_y);
            alpha(3,4) = phi_34(1 - csi_x, 2 - csi_y);
            alpha(4,1) = phi_44(2 - csi_x, 1 + csi_y);
            alpha(4,2) = phi_43(2 - csi_x, csi_y);
            alpha(4,3) = phi_43(2 - csi_x, 1 - csi_y);
            alpha(4,4) = phi_44(2 - csi_x, 2 - csi_y);
            
            % Assign to design matrix
            for k = -1:2
                for h = -1:2
                    if (i_x + k >= 0 && i_x + k <= length(knots_x) - 1 && ...
                        i_y + h >= 0 && i_y + h <= length(knots_y) - 1)
                        % Column index (column-major, matching order(i_x, i_y) = i_y + i_x * yNum)
                        col = (i_y + h) + (i_x + k) * length(knots_y) + 1;
                        A_spline_obs(n, col) = alpha(k + 2, h + 2);
                    end
                end
            end
        end
    end

% Compute Spline Parameter Covariance
N_spline = A_spline_obs' * A_spline_obs;
lambda_mat = lambda_spl * eye(size(N_spline));   % regularization if used
inv_N_spline = inv(N_spline + lambda_mat);

% 5.3.2) Build A_spline_est (prediction Grid)
x_est_vec = xy_est(:,1);
y_est_vec = xy_est(:,2);
A_spline_est = zeros(length(x_est_vec), n_params_spl);

    for n = 1:length(x_est_vec)
        x = x_est_vec(n);
        y = y_est_vec(n);
        
        % Find knot interval for x
        i_x = find(x >= knots_x(1:end-1) & x < knots_x(2:end), 1);
        if isempty(i_x)
            if x <= knots_x(1), i_x = 1; elseif x >= knots_x(end), i_x = length(knots_x) - 1; end
        else
            i_x = i_x - 1; 
        end
        
        % Find knot interval for y
        i_y = find(y <= knots_y(1:end-1) & y > knots_y(2:end), 1);
        if isempty(i_y)
            if y >= knots_y(1), i_y = 1; elseif y <= knots_y(end), i_y = length(knots_y) - 1; end
        else
            i_y = i_y - 1; 
        end
        
        if (i_x >= -2 && i_x <= length(knots_x) && i_y >= -2 && i_y <= length(knots_y))
            csi_x = (x - knots_x(i_x + 1)) / deltaX;
            csi_y = (y - knots_y(i_y + 1)) / deltaY;
            
            alpha = zeros(4, 4);
            alpha(1,1) = phi_44(1 + csi_x, 1 + csi_y); alpha(1,2) = phi_43(1 + csi_x, csi_y); alpha(1,3) = phi_43(1 + csi_x, 1 - csi_y); alpha(1,4) = phi_44(1 + csi_x, 2 - csi_y);
            alpha(2,1) = phi_34(csi_x, 1 + csi_y);     alpha(2,2) = phi_33(csi_x, csi_y);     alpha(2,3) = phi_33(csi_x, 1 - csi_y);     alpha(2,4) = phi_34(csi_x, 2 - csi_y);
            alpha(3,1) = phi_34(1 - csi_x, 1 + csi_y); alpha(3,2) = phi_33(1 - csi_x, csi_y); alpha(3,3) = phi_33(1 - csi_x, 1 - csi_y); alpha(3,4) = phi_34(1 - csi_x, 2 - csi_y);
            alpha(4,1) = phi_44(2 - csi_x, 1 + csi_y); alpha(4,2) = phi_43(2 - csi_x, csi_y); alpha(4,3) = phi_43(2 - csi_x, 1 - csi_y); alpha(4,4) = phi_44(2 - csi_x, 2 - csi_y);
            
            for k = -1:2
                for h = -1:2
                    if (i_x + k >= 0 && i_x + k <= length(knots_x) - 1 && i_y + h >= 0 && i_y + h <= length(knots_y) - 1)
                        col = (i_y + h) + (i_x + k) * length(knots_y) + 1;
                        A_spline_est(n, col) = alpha(k + 2, h + 2);
                    end
                end
            end
        end
    end

% Propagate spline variance
% var = s0^2 * diag(A * inv(N) * A')
var_spline_grid = s02_spl * sum((A_spline_est * inv_N_spline) .* A_spline_est, 2);


% 5.4) Total Variance and Reshape
% (assuming independence between Poly and Spline errors)
var_total_col =  var_spline_grid;

% Reshape to [n_s, n_t]
final_std = reshape(sqrt(var_total_col), length(t_full), length(s_full))';

disp('Variance propagation completed.');




%% 6) OUTPUT VARIABLES
disp('======== Step 6 ========');
disp('Output variables preparation started...');

% 6.1) Define latitude and longitude of output grid
grid_x = centerline_data.xy_centerline(:,1);
grid_y = centerline_data.xy_centerline(:,2);
[lat_centerline, lon_centerline] = utm2deg(grid_x, grid_y, repmat(utmZone, length(grid_x), 1));
lonlatIN_AOI_DET1D = [lon_centerline, lat_centerline];

% Lat/Lon coordinaates for PS plot
[lat_ps, lon_ps] = utm2deg(xyIN_AOI(:, 1), xyIN_AOI(:, 2), repmat(utmZone, size(xyIN_AOI, 1), 1));


% 6.2) Export final_signal_out for estimation grid points
final_signal_out = final_signal_shift;
final_std_out = final_std;

% verify sizes match
assert(size(final_signal_out, 1) == size(lonlatIN_AOI_DET1D, 1), ...
    'Mismatch: final_signal_out has %d rows, lonlatIN_AOI_DET2D has %d rows.', ...
    size(final_signal_out, 1), size(lonlatIN_AOI_DET1D, 1));
assert(size(final_signal_out, 2) == length(t_full), ...
    'Mismatch: final_signal_out has %d columns, expected %d times.', ...
    size(final_signal_out, 2), length(t_full));


% 6.3) Define the final time
dates_full = t_dateIN(1) + t_full;


% 6.4) Re-sort export variable final_signal_orig to match the input one
[~, unsort_idx] = sort(idx_SAR);
final_signal_orig = final_signal_orig(unsort_idx, :);


% OUTPUT:
% lonlatIN_AOI_DET2D = lonlatIN_AOI_DET1D;
% dates_full = dates_full;
% t_full = t_full;
% final_signal_orig = final_signal_orig;
% final_signal = final_signal_out;
% final_std = final_std_out;

disp('Output variables preparation completed.');




%% 7) PLOTS
disp('======== Step 7 ========');
disp('Plots creation started...');

% 7.1) GIF with basemap
figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
v_min = min(round(prctile(final_signal_out(:), 5), 0), -5);
v_max = max(round(prctile(final_signal_out(:), 95), 0), 5);
h = waitbar(0, 'Plotting epochs...');
for t = 1:length(t_full)
    waitbar(t/length(t_full), h, sprintf('Plotting epoch %d/%d', t, length(t_full)));

    displ_at_t = final_signal_shift(:, t);
    displ_at_t_shp = displ_at_t(:);

    geobasemap satellite;
    if t == 1
        pause(20)
    end
    hold on;
    geoscatter(lat_centerline, lon_centerline, markerSize_DET1D, displ_at_t_shp, 'filled');
    geoscatter(lat_ps, lon_ps, markerSize_DET1D/4, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none');

    colormap(jet); clim([v_min, v_max]);
    c = colorbar; c.Label.String = 'LOS Displacement [mm]'; c.Label.FontSize = 15;
    title(sprintf('LOS Displacement on %s (1D)', datestr(dates_full(t))), 'FontSize', 18);
    hold off;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if t == 1
        imwrite(imind, cm, fullfile(figsDir, 'STdet_displ1D.gif'), 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, fullfile(figsDir, 'STdet_displ1D.gif'), 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
    clf;
end
close(h); close;


% 7.2) Heatmap
figure('Visible', 'off', 'Position', [0, 0, 1200, 800]);
imagesc(dates_full, s_full, final_signal_out, 'AlphaData', ~isnan(final_signal_out));
colormap(jet);
ylim([min(s_full), max(s_full)])
ylabel('Distance [m]');
xlabel('Time [date]');
clim([v_min v_max]);
c = colorbar;
c.Label.String = 'LOS displacement'; 
c.Label.FontSize = 15; 
title(c, '[mm]', 'FontSize', 15);
set(gca, 'FontSize', 15)
filename = strcat(figsDir, filesep, 'STdet_displ1D.png');
print(gcf, filename, '-dpng', '-r300')
close

disp('GIF creation completed.');




%%             -------- END OF SCRIPT -------
disp('==================================================');
disp('--- 1D Deterministic-based modelling COMPLETED ---');
disp('==================================================');