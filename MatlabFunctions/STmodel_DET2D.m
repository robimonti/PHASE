function [lonlatIN_AOI_DET2D, dates_full, t_full, final_signal_orig, final_signal_out, final_std_out] = ...
    STmodel_DET2D(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, x_grid, y_grid, figsDir, minMonths, ...
    gS_input_path, gS_output_path, gS_job_path, detectedOS, varNoise_method, varNoise_manual, num_spl_method, ...
    spline_method, num_spl_row_manual, num_spl_col_manual, num_spl_time_manual, lambda_method, lambda_manual, ...
    utmZone, xyIN_AOI, xyAOI, step_t, detrend_method, use_inclined_means, poly_degree, markerSize_DET2D)

% STmodel_DET1D Performs spatio-temporal deterministic modelling for displacement data
%
% Inputs:
%   displIN_AOI           - n x m matrix of displacement data (n points, m epochs)
%   PSidIN_AOI            - n x 1 vector of point IDs
%   t_dateIN              - m x 1 vector of absolute dates (datetime)
%   t_relIN               - m x 1 vector of relative times (days)
%   x_grid                - matrix containing x-coords at all query points
%   y_grid                - matrix containing y-coords at all query points
%   figsDir               - Directory path for saving figures
%   file_out              - Minimum period length for splines interpolation
%   gS_input_path         - path to input folder for geoSplinter
%   gS_output_path        - path to output folder for geoSplinter
%   gS_job_path           - path to job folder for geoSplinter
%   detectedOS            - detected operative system
%   varNoise_method       - a-priori noise variance method 'auto'/'manual'
%   varNoise_manual       - a-priori noise variance value for manual selection (otherwise NaN)
%   num_spl_method        - method for splines number selection 'auto'/'manual'
%   spline_method         - method for 'auto' splines selection 'variance'/'MDL'/'F-test'/'chi2_test'
%   num_spl_row_manual    - manual number of splines for rows (otherwise NaN)
%   num_spl_col_manual    - manual number of splines for columns (otherwise NaN)
%   num_spl_time_manual   - manual number of splines for time (otherwise NaN)
%   lambda_method         - method for lambda estimation 'auto'/'manual'
%   lambda_manual         - value of lambda for 'manual' case' (otherwise NaN)
%   utmZone               - UTM zone for reprojection of coordinates
%   xyIN_AOI              - coordinates of PS inside the AOI
%   xyAOI                 - coordiantes of the AOI limits
%   step_t                - time step for the estimation
%   detrend_method        - option to select the processing method: residual atmospheric errors (residualAtm) or clean data (cleanObs)
%   use_inclined_means    - flag for using flat or inclines plane for residual atmospheric removal
%   poly_degree           - maximum polynomial degree
%   markerSize_DET2D      - marker size for the GIF plot
% 
%
% Output:
%   lonlatIN_AOI_DET2D    - lon/lat coordinates of estimated centerline points
%   dates_full            - dates prediction
%   t_full                - relative time prediction
%   final_signal_orig     - modelled signal in observation coordinates
%   final_signal_out      - modelled signal in query coordinates
%   final_std_out         - uncertainty


disp('--- 2D Deterministic-based modelling started... ---');




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
            gS_exec = fullfile(gS_dir, 'geoSplinter_analysis');
            gS_job_file = fullfile('.', gS_job_path, [file_out, '.job']);
            if isunix
                job_execution = sprintf('%s < %s', gS_exec, gS_job_file);
            else
                temp_bat = [tempname() '.bat'];
                fid = fopen(temp_bat, 'w');
                fprintf(fid, '@echo off\r\n"%s" < "%s"\r\n', gS_exec, gS_job_file);
                fclose(fid);
                job_execution = ['"' temp_bat '"'];
            end
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
        Z_lim = 1.96;
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
displIN_AOI_sorted = displIN_AOI;
data_reg_spl = data_reg_spl';
data_reg_spl_sorted = data_reg_spl;
dataRawCleanT_sorted = data_raw_cleanT;
t_reg_sorted = t_reg_spl;

% Compute residuals of observations - splines interpolation
res_spl = displIN_AOI_sorted - data_reg_spl_sorted;

% Define dimensions
np = size(xyIN_AOI, 1);
nt = length(t_relIN);

% Reshape cleaned observations to column vector (ordered by time, then point)
dataRawCleanT_sorted_col = reshape(dataRawCleanT_sorted', [], 1);
dataRawCleanT_sorted_col(isnan(dataRawCleanT_sorted_col)) = [];

% Create the full grid of observed positions (x, y) and times
grid_t0 = repmat(t_relIN, np, 1);         % [np x nt]
grid_x0 = repmat(xyIN_AOI(:, 1), 1, nt);  % [np x nt]
grid_y0 = repmat(xyIN_AOI(:, 2), 1, nt);  % [np x nt]

% Filter out missing observations
mask_obs = isnan(dataRawCleanT_sorted);  % [np x nt]
grid_t_flt = grid_t0;
grid_t_flt(mask_obs) = NaN;
grid_x_flt = grid_x0;
grid_x_flt(mask_obs) = NaN;
grid_y_flt = grid_y0;
grid_y_flt(mask_obs) = NaN;

% Make them column vectors (same ordering: time-major, then point)
mask_obs_col = reshape(mask_obs', [], 1);  % [nt * np x 1]
grid_t_flt_col = reshape(grid_t_flt', [], 1);
grid_t_flt_col(isnan(grid_t_flt_col)) = [];
grid_x_flt_col = reshape(grid_x_flt', [], 1);
grid_x_flt_col(isnan(grid_x_flt_col)) = [];
grid_y_flt_col = reshape(grid_y_flt', [], 1);
grid_y_flt_col(isnan(grid_y_flt_col)) = [];

% Do the same for the full grid without removed nans
t_full_0 = reshape(grid_t0', [], 1);  % [nt * np x 1]
x_full_0 = reshape(grid_x0', [], 1);
y_full_0 = reshape(grid_y0', [], 1);

residuals_grid_full = dataRawCleanT_sorted;

disp('Variables preparation completed.');




%% 3) LEAST SQUARES MODEL REDUCTION
disp('======== Step 3 ========');
disp('Least Squares model reduction started...');

% observations
y0_lsa = reshape(residuals_grid_full', [], 1);
y0_lsa(isnan(y0_lsa)) = [];

% Prepare the environment
% centering of coordinates to fix ill-conditioning
x0 = mean(grid_x_flt_col);
y0 = mean(grid_y_flt_col);
t0 = mean(grid_t_flt_col);
    
    % centered versions for fitting (observed points)
    x_centered = grid_x_flt_col - x0;
    y_centered = grid_y_flt_col - y0;
    t_centered = grid_t_flt_col - t0;
    % centered versions for reconstruction on the full regular grid
    x_full_centered = x_full_0 - x0;
    y_full_centered = y_full_0 - y0;
    t_full_centered = t_full_0 - t0;


% 3.1) Covariance matrix Q_ll (same for trend and spatial means)
variances = ones(length(y0_lsa), 1);
idx_master = (grid_t_flt_col == master_t_rel);
if any(idx_master)
    variances(idx_master) = 0.1;
end
Q_ll = spdiags(variances, 0, length(y0_lsa), length(y0_lsa));


% 3. Design matrix construction
% ------------------------------------------------------------
alpha_sig = 0.05;
n_epochs = length(t_relIN);

switch detrend_method

    case 'residualAtm'
    % 3.2.1) Residual errors in the data (Orbits + Atmosphere)

    % Global orbital trend (static plane a*x + b*y)
    A_trend = [x_centered, y_centered];
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
        
        % 3 parameters per epoch: [Bias, SlopeX, SlopeY]
        % Total atmos params = 3 * n_epochs
        
        % For each observation (i, t) -> columns corresponding to epoch t.
        % The row is: [ ... 0 0 1 x_i y_i 0 0 ... ] at the columns for epoch t.
        
        % Map each valid observation to its epoch index (1 to n_epochs)
        [~, obs_epoch_idx] = ismember(grid_t_flt_col, t_relIN);
        
        % Base indices for the 3 parameters of each epoch
        % Cols for epoch k are: (k-1)*3 + [1, 2, 3]
        col_bias  = (obs_epoch_idx - 1) * 3 + 1;
        col_slopeX = (obs_epoch_idx - 1) * 3 + 2;
        col_slopeY = (obs_epoch_idx - 1) * 3 + 3;
        
        n_total_obs = length(grid_t_flt_col);
        n_atmos_params = 3 * n_epochs;
        
        % Create Sparse Matrix (Memory efficient for large N)
        % Rows: 1 to N, Cols: computed above, Values: 1, x, y
        I = repmat((1:n_total_obs)', 3, 1);
        J = [col_bias; col_slopeX; col_slopeY];
        V = [ones(n_total_obs, 1); x_centered; y_centered];
        
        A_atmos = sparse(I, J, V, n_total_obs, n_atmos_params);
        
        % Convert to full if N is small enough to speed up solvers, 
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
        rhs = W_A' * y0_lsa;
        
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
        v = y0_lsa - A_curr * x_curr;
        
        % Weighted Sum of Squared Residuals (WSSR) (v' * inv(Q) * v)
        W_v = Q_ll \ v;
        WSSR = v' * W_v;
        
        % Degrees of Freedom
        n_params_curr = sum(mask_params);
        dof = length(y0_lsa) - n_params_curr;
        
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
        A_red_trend = [x_full_centered, y_full_centered];
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
            % x_atmos_all is [3*n_epochs x 1]
            % Reshape to [3 x n_epochs] -> Row 1: Bias, Row 2: SlopeX, Row 3: SlopeY
            atmos_params = reshape(x_atmos_all, 3, n_epochs);
            
            biases  = atmos_params(1, :); % [1 x n_epochs]
            slopesX = atmos_params(2, :);
            slopesY = atmos_params(3, :);
            
            % Construct the full grid correction
            % Atmos(t,x,y) = Bias(t) + Sx(t)*x + Sy(t)*y
            
            % Grid X and Y (centered)
            X_grid = reshape(x_full_centered, length(t_relIN), [])'; % [n_PS x n_epochs]
            Y_grid = reshape(y_full_centered, length(t_relIN), [])';
            
            % Expand parameters to grid size
            B_grid = repmat(biases, size(residuals_grid_full, 1), 1);
            Sx_grid = repmat(slopesX, size(residuals_grid_full, 1), 1);
            Sy_grid = repmat(slopesY, size(residuals_grid_full, 1), 1);
            
            spatialMeans_grid = B_grid + (Sx_grid .* X_grid) + (Sy_grid .* Y_grid);
        end
    
    end

end


% 3.5) Polynomial surface fitting (empirical)
% If 'residualAtm' was run, the polynomial is fitted on the residuals of that model.
% If 'cleanObs' was selected, it is fitted on the raw y0_lsa.
if exist('x_est_all', 'var') && exist('A_full', 'var')
    % remove the previously estimated Orbit + Atmosphere
    y_for_poly = y0_lsa - A_full * x_est_all;
else
    % no previous subtraction
    y_for_poly = y0_lsa;
end

% 3.5.1) Build design matrix A_poly
% (using the vectors corresponding to current VALID observations)
t_n = t_centered;
x_n = x_centered;
y_n = y_centered;

% Start with intercept
A_poly = ones(length(y0_lsa), 1);
param_names = {'Bias'};

% Degree 1
if poly_degree >= 1
    A_poly = [A_poly, x_n, y_n, t_n];
    param_names = [param_names, {'x', 'y', 't'}];
end

% Degree 2
if poly_degree >= 2
    A_poly = [A_poly, x_n.^2, y_n.^2, t_n.^2, ...
                  x_n.*y_n, x_n.*t_n, y_n.*t_n];
    param_names = [param_names, {'x^2', 'y^2', 't^2', 'xy', 'xt', 'yt'}];
end

% Degree 3
if poly_degree >= 3
    A_poly = [A_poly, x_n.^3, y_n.^3, t_n.^3, ...
                  x_n.^2.*y_n, x_n.^2.*t_n, ...
                  y_n.^2.*x_n, y_n.^2.*t_n, ...
                  t_n.^2.*x_n, t_n.^2.*y_n, ...
                  x_n.*y_n.*t_n];
    param_names = [param_names, {'x^3', 'y^3', 't^3', 'x^2y', 'x^2t', ...
                        'y^2x', 'y^2t', 't^2x', 't^2y', 'xyt'}];
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
x_f = x_full_centered;   % [np*nt x 1]
y_f = y_full_centered;   % [np*nt x 1]
t_f = t_full_centered;   % [np*nt x 1]

% Build A_poly_full
A_poly_full = ones(length(x_f), 1);

if poly_degree >= 1
    A_poly_full = [A_poly_full, x_f, y_f, t_f];
end

if poly_degree >= 2
    A_poly_full = [A_poly_full, x_f.^2, y_f.^2, t_f.^2, ...
                       x_f.*y_f, x_f.*t_f, y_f.*t_f];
end

if poly_degree >= 3
    A_poly_full = [A_poly_full, x_f.^3, y_f.^3, t_f.^3, ...
                       x_f.^2.*y_f, x_f.^2.*t_f, ...
                       y_f.^2.*x_f, y_f.^2.*t_f, ...
                       t_f.^2.*x_f, t_f.^2.*y_f, ...
                       x_f.*y_f.*t_f];
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
x_coords = reshape(repmat(xyIN_AOI(:,1)', nt, 1), [], 1);  % [np*nt x 1], time-major
y_coords = reshape(repmat(xyIN_AOI(:,2)', nt, 1), [], 1);
t_coords = repmat(t_relIN, np, 1);                         % [np x nt], then flatten to [np*nt x 1]
t_coords = reshape(t_coords', [], 1);
obs = reshape(data_final_clean', [], 1);                   % [np*nt x 1]

% remove NaNs for fitting
valid_idx = ~isnan(obs);
x_coords = x_coords(valid_idx);
y_coords = y_coords(valid_idx);
t_coords = t_coords(valid_idx);
obs = obs(valid_idx);
n_obs = length(obs);

    % 4.2.1) Ranges for normalization
    x_min_val = min(xyIN_AOI(:,1)); 
    x_max_val = max(xyIN_AOI(:,1));
    y_min_val = min(xyIN_AOI(:,2)); 
    y_max_val = max(xyIN_AOI(:,2));
    t_min_val = min(t_relIN);       
    t_max_val = max(t_relIN);
    
    x_range = x_max_val - x_min_val;
    y_range = y_max_val - y_min_val;
    t_range = t_max_val - t_min_val;
    
    % Define scaling factors (Ranges)
    % Use the maximum spatial dimension to preserve X/Y aspect ratio (isotropy)
    range_s = max(x_range, y_range); 
    range_t = t_range;

    % scale coordinates using correlation ranges for isotropic dimensionless space
    x_scaled = x_coords / range_s;
    y_scaled = y_coords / range_s;
    t_scaled = t_coords / range_t;
    
    % points and values
    points = [x_scaled, y_scaled, t_scaled];
    values = obs;
    
    % find nearest neighbors using knnsearch for efficiency
    [idx, min_dists] = knnsearch(points, points, 'K', 2);
    closest_idxs = idx(:,2);
    min_dists = min_dists(:,2);
    
    % valid mask (avoid zero dist)
    valid_mask = min_dists > 1e-10;
    values_i = values(valid_mask);
    values_closest = values(closest_idxs(valid_mask));
    points_valid = points(valid_mask,:);
    points_closest = points(closest_idxs(valid_mask),:);
    
    % first derivatives and midpoints
    derivatives = (values_closest - values_i) ./ min_dists(valid_mask);
    midpoints = (points_valid + points_closest) / 2;
    
    % unique midpoints (mimic 2D code, though rare duplicates in 3D)
    [unique_midpoints, unique_idx] = unique(midpoints, 'rows', 'stable');
    derivatives = derivatives(unique_idx);
    
    % nearest for unique midpoints
    [idx_mid, min_dists_mid] = knnsearch(unique_midpoints, unique_midpoints, 'K', 2);
    closest_idxs_mid = idx_mid(:,2);
    min_dists_mid = min_dists_mid(:,2);
    
    % valid mask for mid
    valid_mask_mid = min_dists_mid > 1e-10;
    derivatives_i = derivatives(valid_mask_mid);
    derivatives_closest = derivatives(closest_idxs_mid(valid_mask_mid));
    
    % second derivatives
    new_derivatives = (derivatives_closest - derivatives_i) ./ min_dists_mid(valid_mask_mid);
    
    % sigma2n (fallback if empty/invalid)
    if isempty(new_derivatives) || all(~isfinite(new_derivatives))
        sigma2n = 1e-3;
    else
        sigma2n = mean(new_derivatives.^2, 'omitnan');
        if ~isfinite(sigma2n) || sigma2n <= 0
            sigma2n = 1e-3;
        end
    end


% 4.3) Processing domains
% Initial estimate for num_row and num_col based on correlation ranges
base_num_x = max(2, round(x_range / range_s));
base_num_y = max(2, round(y_range / range_s));
base_num_t = max(2, round(t_range / range_t));

% Define min and max number of splines to test in each direction
min_num_x = max(2, round(base_num_x * 0.5));            % 50% of base
max_num_x = min(round(base_num_x * 2), round(np / 4));  % 200% or 1/4 points
min_num_y = max(2, round(base_num_y * 0.5));
max_num_y = min(round(base_num_y * 2), round(np / 4));
min_num_t = max(2, round(base_num_t * 0.5));
max_num_t = min(round(base_num_t * 2), round(nt / 4));

% Test values (limit combinations)
x_test = unique(round(linspace(min_num_x, max_num_x, 5)));
y_test = unique(round(linspace(min_num_y, max_num_y, 5)));
t_test = unique(round(linspace(min_num_t, max_num_t, 5)));


% 4.4) Initialize variables for candidate solution
switch num_spl_method
    case 'manual'
        num_x = num_spl_row_manual;
        num_y = num_spl_col_manual;
        num_t = num_spl_time_manual;
        combinations = [num_x, num_y, num_t];
    case 'auto'
        [X, Y, T] = meshgrid(x_test, y_test, t_test);
        combinations = [X(:), Y(:), T(:)];
end

% Pre-allocate for candidates
num_candidates = size(combinations, 1);
obs_est_tmp = NaN(n_obs, num_candidates);
var_resSpl = NaN(num_candidates, 1);
s02_tmp = NaN(num_candidates, 1);
MDL_tmp = NaN(num_candidates, 1);
combs = NaN(num_candidates, 3);

% Basis functions
phi_3 = @(csi) ((2 - csi).^3 - 4 * (1 - csi).^3) / 6;
phi_4 = @(csi) (2 - csi).^3 / 6;

% Define 3D tensor-product basis functions
phi_333 = @(csi_x, csi_y, csi_t) phi_3(csi_x) .* phi_3(csi_y) .* phi_3(csi_t);
phi_334 = @(csi_x, csi_y, csi_t) phi_3(csi_x) .* phi_3(csi_y) .* phi_4(csi_t);
phi_343 = @(csi_x, csi_y, csi_t) phi_3(csi_x) .* phi_4(csi_y) .* phi_3(csi_t);
phi_344 = @(csi_x, csi_y, csi_t) phi_3(csi_x) .* phi_4(csi_y) .* phi_4(csi_t);
phi_433 = @(csi_x, csi_y, csi_t) phi_4(csi_x) .* phi_3(csi_y) .* phi_3(csi_t);
phi_434 = @(csi_x, csi_y, csi_t) phi_4(csi_x) .* phi_3(csi_y) .* phi_4(csi_t);
phi_443 = @(csi_x, csi_y, csi_t) phi_4(csi_x) .* phi_4(csi_y) .* phi_3(csi_t);
phi_444 = @(csi_x, csi_y, csi_t) phi_4(csi_x) .* phi_4(csi_y) .* phi_4(csi_t);


% 4.5) Iterate over candidates
for cc = 1:num_candidates

    % 4.5.1) Define dimensions based on current candidate
    num_x = combinations(cc, 1);
    num_y = combinations(cc, 2);
    num_t = combinations(cc, 3);
    
    % Control points count
    n_ctrl_x = num_x + 2; 
    n_ctrl_y = num_y + 2; 
    n_ctrl_t = num_t + 2;
    n_params = n_ctrl_x * n_ctrl_y * n_ctrl_t;
    
    if n_obs <= n_params
        warning('n_obs (%d) <= n_params (%d). Skipping.', n_obs, n_params);
        continue;
    end
    
    % 4.5.2) Define knots strictly on the AOI
    deltaX = x_range / num_x;
    deltaY = y_range / num_y;
    deltaT = t_range / num_t;
    
    knots_x = min(xyIN_AOI(:,1)) + (0:num_x) * deltaX;
    knots_y = min(xyIN_AOI(:,2)) + (0:num_y) * deltaY;
    knots_t = min(t_relIN) + (0:num_t) * deltaT;
    
    % 4.5.3) Build design matrix A
    A = sparse(n_obs, n_params);
    
    % Pre-compute inverse deltas for speed
    inv_dX = 1/deltaX;
    inv_dY = 1/deltaY;
    inv_dT = 1/deltaT;

    for n = 1:n_obs
        x = x_coords(n);
        y = y_coords(n);
        t = t_coords(n);
        
        % --- X Interval ---
        i_x = floor((x - knots_x(1)) * inv_dX);
        % handle boundary case (right edge)
        if i_x >= num_x, i_x = num_x - 1; end
        if i_x < 0, i_x = 0; end % Safety
        
        % --- Y Interval ---
        i_y = floor((y - knots_y(1)) * inv_dY);
        if i_y >= num_y, i_y = num_y - 1; end
        if i_y < 0, i_y = 0; end
        
        % --- T Interval ---
        i_t = floor((t - knots_t(1)) * inv_dT);
        if i_t >= num_t, i_t = num_t - 1; end
        if i_t < 0, i_t = 0; end

        % Compute Local Coordinates (csi) in [0,1]
        % note: i_x+1 is the knot at the start of the interval
        csi_x = (x - knots_x(i_x + 1)) * inv_dX;
        csi_y = (y - knots_y(i_y + 1)) * inv_dY; 
        csi_t = (t - knots_t(i_t + 1)) * inv_dT;
        
        % Compute alpha coefficients (4x4x4)
        alpha = zeros(4, 4, 4);
        
        % -- slice 1: l = -1 (corresponds to alpha(:,:,1)) --
        % k = -1 (alpha(1,:,1))
        alpha(1,1,1) = phi_444(1+csi_x, 1+csi_y, 1+csi_t);
        alpha(1,2,1) = phi_434(1+csi_x,   csi_y, 1+csi_t);
        alpha(1,3,1) = phi_434(1+csi_x, 1-csi_y, 1+csi_t);
        alpha(1,4,1) = phi_444(1+csi_x, 2-csi_y, 1+csi_t);
        % k = 0 (alpha(2,:,1))
        alpha(2,1,1) = phi_344(  csi_x, 1+csi_y, 1+csi_t);
        alpha(2,2,1) = phi_334(  csi_x,   csi_y, 1+csi_t);
        alpha(2,3,1) = phi_334(  csi_x, 1-csi_y, 1+csi_t);
        alpha(2,4,1) = phi_344(  csi_x, 2-csi_y, 1+csi_t);
        % k = 1 (alpha(3,:,1))
        alpha(3,1,1) = phi_344(1-csi_x, 1+csi_y, 1+csi_t);
        alpha(3,2,1) = phi_334(1-csi_x,   csi_y, 1+csi_t);
        alpha(3,3,1) = phi_334(1-csi_x, 1-csi_y, 1+csi_t);
        alpha(3,4,1) = phi_344(1-csi_x, 2-csi_y, 1+csi_t);
        % k = 2 (alpha(4,:,1))
        alpha(4,1,1) = phi_444(2-csi_x, 1+csi_y, 1+csi_t);
        alpha(4,2,1) = phi_434(2-csi_x,   csi_y, 1+csi_t);
        alpha(4,3,1) = phi_434(2-csi_x, 1-csi_y, 1+csi_t);
        alpha(4,4,1) = phi_444(2-csi_x, 2-csi_y, 1+csi_t);
        
        % -- slice 2: l = 0 (corresponds to alpha(:,:,2)) --
        % k = -1 (alpha(1,:,2))
        alpha(1,1,2) = phi_443(1+csi_x, 1+csi_y,   csi_t);
        alpha(1,2,2) = phi_433(1+csi_x,   csi_y,   csi_t);
        alpha(1,3,2) = phi_433(1+csi_x, 1-csi_y,   csi_t);
        alpha(1,4,2) = phi_443(1+csi_x, 2-csi_y,   csi_t);
        % k = 0 (alpha(2,:,2))
        alpha(2,1,2) = phi_343(  csi_x, 1+csi_y,   csi_t);
        alpha(2,2,2) = phi_333(  csi_x,   csi_y,   csi_t);
        alpha(2,3,2) = phi_333(  csi_x, 1-csi_y,   csi_t);
        alpha(2,4,2) = phi_343(  csi_x, 2-csi_y,   csi_t);
        % k = 1 (alpha(3,:,2))
        alpha(3,1,2) = phi_343(1-csi_x, 1+csi_y,   csi_t);
        alpha(3,2,2) = phi_333(1-csi_x,   csi_y,   csi_t);
        alpha(3,3,2) = phi_333(1-csi_x, 1-csi_y,   csi_t);
        alpha(3,4,2) = phi_343(1-csi_x, 2-csi_y,   csi_t);
        % k = 2 (alpha(4,:,2))
        alpha(4,1,2) = phi_443(2-csi_x, 1+csi_y,   csi_t);
        alpha(4,2,2) = phi_433(2-csi_x,   csi_y,   csi_t);
        alpha(4,3,2) = phi_433(2-csi_x, 1-csi_y,   csi_t);
        alpha(4,4,2) = phi_443(2-csi_x, 2-csi_y,   csi_t);
        
        % -- slice 3: l = 1 (corresponds to alpha(:,:,3)) --
        % k = -1 (alpha(1,:,3))
        alpha(1,1,3) = phi_443(1+csi_x, 1+csi_y, 1-csi_t);
        alpha(1,2,3) = phi_433(1+csi_x,   csi_y, 1-csi_t);
        alpha(1,3,3) = phi_433(1+csi_x, 1-csi_y, 1-csi_t);
        alpha(1,4,3) = phi_443(1+csi_x, 2-csi_y, 1-csi_t);
        % k = 0 (alpha(2,:,3))
        alpha(2,1,3) = phi_343(  csi_x, 1+csi_y, 1-csi_t);
        alpha(2,2,3) = phi_333(  csi_x,   csi_y, 1-csi_t);
        alpha(2,3,3) = phi_333(  csi_x, 1-csi_y, 1-csi_t);
        alpha(2,4,3) = phi_343(  csi_x, 2-csi_y, 1-csi_t);
        % k = 1 (alpha(3,:,3))
        alpha(3,1,3) = phi_343(1-csi_x, 1+csi_y, 1-csi_t);
        alpha(3,2,3) = phi_333(1-csi_x,   csi_y, 1-csi_t);
        alpha(3,3,3) = phi_333(1-csi_x, 1-csi_y, 1-csi_t);
        alpha(3,4,3) = phi_343(1-csi_x, 2-csi_y, 1-csi_t);
        % k = 2 (alpha(4,:,3))
        alpha(4,1,3) = phi_443(2-csi_x, 1+csi_y, 1-csi_t);
        alpha(4,2,3) = phi_433(2-csi_x,   csi_y, 1-csi_t);
        alpha(4,3,3) = phi_433(2-csi_x, 1-csi_y, 1-csi_t);
        alpha(4,4,3) = phi_443(2-csi_x, 2-csi_y, 1-csi_t);
        
        % -- slice 4: l = 2 (corresponds to alpha(:,:,4)) --
        % k = -1 (alpha(1,:,4))
        alpha(1,1,4) = phi_444(1+csi_x, 1+csi_y, 2-csi_t);
        alpha(1,2,4) = phi_434(1+csi_x,   csi_y, 2-csi_t);
        alpha(1,3,4) = phi_434(1+csi_x, 1-csi_y, 2-csi_t);
        alpha(1,4,4) = phi_444(1+csi_x, 2-csi_y, 2-csi_t);
        % k = 0 (alpha(2,:,4))
        alpha(2,1,4) = phi_344(  csi_x, 1+csi_y, 2-csi_t);
        alpha(2,2,4) = phi_334(  csi_x,   csi_y, 2-csi_t);
        alpha(2,3,4) = phi_334(  csi_x, 1-csi_y, 2-csi_t);
        alpha(2,4,4) = phi_344(  csi_x, 2-csi_y, 2-csi_t);
        % k = 1 (alpha(3,:,4))
        alpha(3,1,4) = phi_344(1-csi_x, 1+csi_y, 2-csi_t);
        alpha(3,2,4) = phi_334(1-csi_x,   csi_y, 2-csi_t);
        alpha(3,3,4) = phi_334(1-csi_x, 1-csi_y, 2-csi_t);
        alpha(3,4,4) = phi_344(1-csi_x, 2-csi_y, 2-csi_t);
        % k = 2 (alpha(4,:,4))
        alpha(4,1,4) = phi_444(2-csi_x, 1+csi_y, 2-csi_t);
        alpha(4,2,4) = phi_434(2-csi_x,   csi_y, 2-csi_t);
        alpha(4,3,4) = phi_434(2-csi_x, 1-csi_y, 2-csi_t);
        alpha(4,4,4) = phi_444(2-csi_x, 2-csi_y, 2-csi_t);

        % Fill Matrix
        for k = -1:2
            for h = -1:2
                for l = -1:2
                    % Indices relative to control points (0-based internal logic)
                    abs_ix = i_x + k + 1; % Shift +1 because B-splines usually -1..2
                    abs_iy = i_y + h + 1;
                    abs_it = i_t + l + 1;

                    % Check bounds
                    % valid range: 1 to n_ctrl
                    if (abs_ix >= 1 && abs_ix <= n_ctrl_x && ...
                        abs_iy >= 1 && abs_iy <= n_ctrl_y && ...
                        abs_it >= 1 && abs_it <= n_ctrl_t)
                        
                        % Linear Index: (y, x, t) order
                        % col = y + (x-1)*Ny + (t-1)*Ny*Nx
                        col = abs_iy + (abs_ix - 1) * n_ctrl_y + (abs_it - 1) * n_ctrl_y * n_ctrl_x;
                        
                        A(n, col) = alpha(k + 2, h + 2, l + 2);
                    end
                end
            end
        end
    end

    % Regularization parameter
    switch lambda_method
        case 'auto'
            % Effective scaled deltaGrid (geometric mean for 3D)
            deltaX_scaled = deltaX / range_s;
            deltaY_scaled = deltaY / range_s;
            deltaT_scaled = deltaT / range_t;
            deltaGrid_scaled = (deltaX_scaled * deltaY_scaled * deltaT_scaled)^(1/3);
            
            lambda = varNoise / (sigma2n * deltaGrid_scaled^4);
        case 'manual'
            lambda = lambda_manual;
    end

    % Solve least squares with regularization
    N = A' * A + lambda * eye(n_params);
    rhs = A' * obs;
    if cond(N) > 1e10
        coef = pinv(N) * rhs;
    else
        coef = N \ rhs;
    end
    estimated_obs = A * coef;

    % Residuals and metrics
    residuals = obs - estimated_obs;
    SSR = sum(residuals.^2);
    dof = n_obs - n_params;
    s02_est = SSR / dof;

    % Store the computed variables
    obs_est_tmp(:, cc) = estimated_obs;
    s02_tmp(cc) = s02_est;
    var_resSpl(cc) = s02_est;
    MDL_tmp(cc) = log(n_obs) * n_params / n_obs + log(SSR + eps);
    lambda_tmp(cc) = lambda;
    combs(cc, :) = [num_x, num_y, num_t];
    coef_tmp{cc} = coef;

end

% Calculate parameters for all candidates
n_params_all = (combs(:,1) + 2) .* (combs(:,2) + 2) .* (combs(:,3) + 2);

% Sort in ascending order (simplest -> most complex)
[n_params_sorted, sort_idx] = sort(n_params_all);

% Reorder all arrays to match this sorted order
combs       = combs(sort_idx, :);
s02_tmp     = s02_tmp(sort_idx);
MDL_tmp     = MDL_tmp(sort_idx);
lambda_tmp  = lambda_tmp(sort_idx);
coef_tmp    = coef_tmp(sort_idx);
obs_est_tmp = obs_est_tmp(:, sort_idx); 

% Only if they exist/calculated
if exist('var_resSpl', 'var'), var_resSpl = var_resSpl(sort_idx); end

var_idx = NaN; % Initialize


% 4.6) Additional metrics computation
switch num_spl_method
    case 'manual'
        % Manual selection: find the user's combination in the sorted list
        num_x = num_spl_row_manual;
        num_y = num_spl_col_manual;
        num_t = num_spl_time_manual;
        
        var_idx = find(combs(:,1) == num_x & ...
                       combs(:,2) == num_y & ...
                       combs(:,3) == num_t, 1);
                   
        if isempty(var_idx)
            error('Manual combination [%d, %d, %d] not found in the tested grid.', ...
                  num_x, num_y, num_t);
        end

    case 'auto'
        switch spline_method
            case 'MDL'
                % MDL is independent of order, just find the minimum
                [~, var_idx] = min(MDL_tmp);

            case 'variance'
                % Variance threshold method
                var_thrs = prctile(var_resSpl, 80);
                var_spl_diff = var_resSpl - var_thrs;
                var_spl_diff(var_spl_diff < 0) = NaN;
                [~, var_idx] = min(var_spl_diff, [], 'omitnan');
                if isempty(var_idx)
                    var_idx = 1; % default to simplest (first in sorted list)
                end

            case 'F_test'
                % Start assuming the simplest model (index 1) is best
                current_best_idx = 1; 
                
                % Iterate through all more complex models
                for kk = 2:num_candidates
                    % Hypothesis: does model 'kk' improve upon 'current_best_idx'?
                    idx_simple  = current_best_idx;
                    idx_complex = kk;
                    
                    dof_simple  = n_obs - n_params_sorted(idx_simple);
                    dof_complex = n_obs - n_params_sorted(idx_complex);
                    
                    SSR_simple  = s02_tmp(idx_simple) * dof_simple;
                    SSR_complex = s02_tmp(idx_complex) * dof_complex;
                    
                    % F-test Validity Check (Must have degrees of freedom diff)
                    if dof_simple > dof_complex && dof_complex > 0
                        
                        numerator   = (SSR_simple - SSR_complex) / (dof_simple - dof_complex);
                        denominator = s02_tmp(idx_complex);
                        F_stat      = numerator / denominator;
                        
                        % p-value: Prob of observing this F if H0 (no improvement) is true
                        p_val = 1 - fcdf(F_stat, dof_simple - dof_complex, dof_complex);
                        
                        % If p < 0.05, the improvement is statistically significant.
                        % We reject H0 and update the current best model.
                        if p_val < alpha_sig
                            current_best_idx = kk;
                        end
                    end
                end
                
                var_idx = current_best_idx;

            case 'chi2_test'
                % Chi2 Test: checks if residuals are consistent with noise variance
                chi2_pvals = NaN(num_candidates, 1);
                
                for kk = 1:num_candidates
                    dof = n_obs - n_params_sorted(kk);
                    if dof > 0
                        chi2_stat = s02_tmp(kk) * dof / varNoise;
                        chi2_pvals(kk) = 1 - chi2cdf(chi2_stat, dof);
                    end
                end
                
                % Look for first model fitting the noise distribution
                idx = find(chi2_pvals >= alpha_sig/2 & chi2_pvals <= (1 - alpha_sig/2), 1, 'first');
                
                if isempty(idx)
                    var_idx = 1; % default to simplest if none fit perfectly
                else
                    var_idx = idx;
                end
                
            otherwise
                error('Invalid spline_method: %s', spline_method);
        end
end


% 4.7) Final retrieval and knots reconstruction
% Retrieve best parameters from the sorted index
num_x = combs(var_idx, 1);
num_y = combs(var_idx, 2);
num_t = combs(var_idx, 3);

n_ctrl_x = num_x + 2;
n_ctrl_y = num_y + 2;
n_ctrl_t = num_t + 2;
n_params = n_ctrl_x * n_ctrl_y * n_ctrl_t;

PS_data_spline = obs_est_tmp(:, var_idx);
res_spl        = obs - PS_data_spline;
lambda_spl     = lambda_tmp(var_idx);
s02_spl        = s02_tmp(var_idx);
coef_spl       = coef_tmp{var_idx};

% Rebuild knots for the selected combination (Increasing Order)
deltaX = x_range / num_x;
deltaY = y_range / num_y;
deltaT = t_range / num_t;

knots_x = min(xyIN_AOI(:,1)) + (0:num_x) * deltaX;
knots_y = min(xyIN_AOI(:,2)) + (0:num_y) * deltaY; 
knots_t = min(t_relIN)       + (0:num_t) * deltaT;


% 4.8) Splines solution estimation on the full query grid
% define full grid of time
t_full = (min(t_relIN):step_t:max(t_relIN))';

% Check if x_grid/y_grid are matrices (from meshgrid) or vectors
if ~isvector(x_grid)
    x_vec = x_grid(1, :)';
    y_vec = y_grid(:, 1);
    if size(x_grid, 1) > 1 && x_grid(1,1) == x_grid(2,1)
         x_vec = unique(x_grid);
         y_vec = unique(y_grid);
    end
else
    x_vec = unique(x_grid(:));
    y_vec = unique(y_grid(:));
end

% Ensure dimensions match
nx = length(x_vec);
ny = length(y_vec);
nt_full = length(t_full);

% Initialize output [ny x nx x nt_full]
% Convention: Row=Y, Col=X, Page=Time
splines_signal = zeros(ny, nx, nt_full);

% Pre-compute Inverse Deltas (Speed)
inv_dX = 1/deltaX;
inv_dY = 1/deltaY;
inv_dT = 1/deltaT;

% 4.8.1) Loop over time epochs to build A
for ti = 1:nt_full
    curr_t = t_full(ti);
    i_t = floor((curr_t - knots_t(1)) * inv_dT); 
    if i_t >= num_t, i_t = num_t - 1; end; if i_t < 0, i_t = 0; end
    csi_t = (curr_t - knots_t(i_t + 1)) * inv_dT;
    B_t = bspline_basis(i_t, csi_t);
    
    slice_2d = zeros(ny, nx);
    for yi = 1:ny
        curr_y = y_vec(yi);
        i_y = floor((curr_y - knots_y(1)) * inv_dY); 
        if i_y >= num_y, i_y = num_y - 1; end; if i_y < 0, i_y = 0; end
        csi_y = (curr_y - knots_y(i_y + 1)) * inv_dY;
        B_y = bspline_basis(i_y, csi_y);
        
        for xi = 1:nx
            curr_x = x_vec(xi);
            i_x = floor((curr_x - knots_x(1)) * inv_dX); 
            if i_x >= num_x, i_x = num_x - 1; end; if i_x < 0, i_x = 0; end
            csi_x = (curr_x - knots_x(i_x + 1)) * inv_dX;
            B_x = bspline_basis(i_x, csi_x);
            
            val = 0;
            for l = -1:2 
                idx_t = i_t + l + 1; % shift +1 for 1-based indexing
                for h = -1:2 
                    idx_y = i_y + h + 1;
                    for k = -1:2 
                        idx_x = i_x + k + 1;
                        
                        if (idx_x >= 1 && idx_x <= n_ctrl_x && ...
                            idx_y >= 1 && idx_y <= n_ctrl_y && ...
                            idx_t >= 1 && idx_t <= n_ctrl_t)
                            
                            coef_idx = idx_y + (idx_x - 1) * n_ctrl_y + (idx_t - 1) * n_ctrl_y * n_ctrl_x;
                            % Access basis with k+2 (mapping -1->1, 0->2, etc.)
                            val = val + coef_spl(coef_idx) * B_x(k+2) * B_y(h+2) * B_t(l+2);
                        end
                    end
                end
            end
            slice_2d(yi, xi) = val;
        end
    end
    splines_signal(:, :, ti) = slice_2d;
end


% 4.9) Restore the removed polynomial surface
% Initialize polynomial signal container
[n_rows, n_cols, n_times] = size(splines_signal);
n_spatial = n_rows * n_cols;
poly_signal = zeros(n_rows, n_cols, n_times);

% Pre-calculate spatial centering (constant for all times)
% x_grid and y_grid are [n_rows x n_cols]
x_n_spatial = x_grid(:) - x0;
y_n_spatial = y_grid(:) - y0;

% Process slice-by-slice to save memory
for ti = 1:n_times
    
    % Current centered time (scalar for this slice)
    t_current_n = t_full(ti) - t0;
    
    % Expand time to match spatial grid size [n_spatial x 1]
    t_vec = repmat(t_current_n, n_spatial, 1);
    
    % Build Design Matrix for this time slice only
    % Base: intercept
    A_slice = ones(n_spatial, 1);
    
    % Degree 1 (linear)
    if poly_degree >= 1
        A_slice = [A_slice, x_n_spatial, y_n_spatial, t_vec];
    end
    
    % Degree 2 (quadratic)
    if poly_degree >= 2
        A_slice = [A_slice, x_n_spatial.^2, y_n_spatial.^2, t_vec.^2, ...
                   x_n_spatial.*y_n_spatial, x_n_spatial.*t_vec, y_n_spatial.*t_vec];
    end
    
    % Degree 3 (cubic)
    if poly_degree >= 3
        A_slice = [A_slice, x_n_spatial.^3, y_n_spatial.^3, t_vec.^3, ...
                   x_n_spatial.^2.*y_n_spatial, x_n_spatial.^2.*t_vec, ...
                   y_n_spatial.^2.*x_n_spatial, y_n_spatial.^2.*t_vec, ...
                   t_vec.^2.*x_n_spatial, t_vec.^2.*y_n_spatial, ...
                   x_n_spatial.*y_n_spatial.*t_vec];
    end
    
    % Estimate slice and reshape
    % x_poly_coeffs is the the full vector (including zeros)
    z_slice = A_slice * x_poly_coeffs;
    
    % Reshape column vector back to [n_rows, n_cols]
    poly_signal(:,:,ti) = reshape(z_slice, n_rows, n_cols);
end


% 4.10) Get the final signal
final_signal = splines_signal + poly_signal;


% 4.11) Apply reference shift (relative to first epoch, t=0)
% This ensures all displacements start at 0 relative to t_start
final_signal_shift   = final_signal   - final_signal(:,:,1);
splines_signal_shift = splines_signal - splines_signal(:,:,1);
poly_signal_shift    = poly_signal    - poly_signal(:,:,1);


% 4.12) Prediction on observation points
% Define target coordinates (full original grid)
% These variables were created in Step 2
target_t = t_full_0; % [np*nt x 1]
target_x = x_full_0;
target_y = y_full_0;

% Initialize output
n_targets = length(target_t);
signal_spl_full = zeros(n_targets, 1);

% Pre-compute inverse deltas for speed
inv_dX = 1/deltaX;
inv_dY = 1/deltaY;
inv_dT = 1/deltaT;

% Direct evaluation loop
% Iteration over every space-time point in the original data structure.
for n = 1:n_targets
    x_val = target_x(n); y_val = target_y(n); t_val = target_t(n);
    
    i_x = floor((x_val - knots_x(1)) * inv_dX); 
    if i_x >= num_x, i_x = num_x - 1; end
    if i_x < 0, i_x = 0; end
    csi_x = (x_val - knots_x(i_x + 1)) * inv_dX;
    i_y = floor((y_val - knots_y(1)) * inv_dY); 
    if i_y >= num_y, i_y = num_y - 1; end
    if i_y < 0, i_y = 0; end
    csi_y = (y_val - knots_y(i_y + 1)) * inv_dY;
    i_t = floor((t_val - knots_t(1)) * inv_dT); 
    if i_t >= num_t, i_t = num_t - 1; end 
    if i_t < 0, i_t = 0; end
    csi_t = (t_val - knots_t(i_t + 1)) * inv_dT;
    
    B_x = bspline_basis(i_x, csi_x); B_y = bspline_basis(i_y, csi_y); B_t = bspline_basis(i_t, csi_t);
    
    val = 0;
    for l = -1:2 
        idx_t = i_t + l + 1;
        for h = -1:2 
            idx_y = i_y + h + 1;
            for k = -1:2 
                idx_x = i_x + k + 1;
                
                if (idx_x >= 1 && idx_x <= n_ctrl_x && ...
                    idx_y >= 1 && idx_y <= n_ctrl_y && ...
                    idx_t >= 1 && idx_t <= n_ctrl_t)
                    
                    coef_idx = idx_y + (idx_x - 1) * n_ctrl_y + (idx_t - 1) * n_ctrl_y * n_ctrl_x;
                    val = val + coef_spl(coef_idx) * B_x(k+2) * B_y(h+2) * B_t(l+2);
                end
            end
        end
    end
    signal_spl_full(n) = val;
end


% 4.12.2) Compute polynomial component at observation points
% centering (as Step 3)
t_target_n = target_t - t0;
x_target_n = target_x - x0;
y_target_n = target_y - y0;

% Build design matrix (A_poly_target)
A_poly_target = ones(length(t_target_n), 1);

if poly_degree >= 1
    A_poly_target = [A_poly_target, x_target_n, y_target_n, t_target_n];
end

if poly_degree >= 2
    A_poly_target = [A_poly_target, x_target_n.^2, y_target_n.^2, t_target_n.^2, ...
                     x_target_n.*y_target_n, x_target_n.*t_target_n, y_target_n.*t_target_n];
end

if poly_degree >= 3
    A_poly_target = [A_poly_target, x_target_n.^3, y_target_n.^3, t_target_n.^3, ...
                     x_target_n.^2.*y_target_n, x_target_n.^2.*t_target_n, ...
                     y_target_n.^2.*x_target_n, y_target_n.^2.*t_target_n, ...
                     t_target_n.^2.*x_target_n, t_target_n.^2.*y_target_n, ...
                     x_target_n.*y_target_n.*t_target_n];
end

% Apply coefficients
signal_poly_full = A_poly_target * x_poly_coeffs;

% 4.12.3) Final summation and reshape
final_signal_orig_col = signal_spl_full + signal_poly_full;
final_signal_orig = reshape(final_signal_orig_col, nt, np)';

disp('Splines interpolation completed.');




%% 5) TOTAL VARIANCE PROPAGATION
disp('======== Step 5 ========');
disp('Variance propagation started...');

% 5.1) Construct design matrices
% Coordinates for observations (centered)
x_obs_n = grid_x_flt_col - x0;
y_obs_n = grid_y_flt_col - y0;
t_obs_n = grid_t_flt_col - t0;

% 5.1.1) Build full A_poly_obs
A_poly_obs = ones(length(x_obs_n), 1);
if poly_degree >= 1
    A_poly_obs = [A_poly_obs, x_obs_n, y_obs_n, t_obs_n];
end
if poly_degree >= 2
    A_poly_obs = [A_poly_obs, x_obs_n.^2, y_obs_n.^2, t_obs_n.^2, ...
                  x_obs_n.*y_obs_n, x_obs_n.*t_obs_n, y_obs_n.*t_obs_n];
end
if poly_degree >= 3
    A_poly_obs = [A_poly_obs, x_obs_n.^3, y_obs_n.^3, t_obs_n.^3, ...
                  x_obs_n.^2.*y_obs_n, x_obs_n.^2.*t_obs_n, ...
                  y_obs_n.^2.*x_obs_n, y_obs_n.^2.*t_obs_n, ...
                  t_obs_n.^2.*x_obs_n, t_obs_n.^2.*y_obs_n, ...
                  x_obs_n.*y_obs_n.*t_obs_n];
end

% Apply the mask from backward elimination
A_poly_obs_red = A_poly_obs(:, mask_poly);

% Compute polynomial parameter covariance matrix (C_xx_poly)
% N = A' * W * A
% C_xx = sigma0^2 * inv(N)
W_A_poly = Q_ll \ A_poly_obs_red;
N_poly   = A_poly_obs_red' * W_A_poly;
if cond(N_poly) > 1e12
    C_xx_poly = sigma0_sq * pinv(N_poly);
else
    C_xx_poly = sigma0_sq * inv(N_poly);
end

% 5.1.2) Define estimation coordinates
x_est_col = repmat(x_grid(:), n_times, 1); 
y_est_col = repmat(y_grid(:), n_times, 1);

% Time coordinates
t_est_col = kron(t_full, ones(n_spatial, 1));

% Center them
x_pred_n = x_est_col - x0;
y_pred_n = y_est_col - y0;
t_pred_n = t_est_col - t0;

% Build A_poly_pred (design matrix on prediction grid)
A_poly_pred = ones(length(x_pred_n), 1);
if poly_degree >= 1
    A_poly_pred = [A_poly_pred, x_pred_n, y_pred_n, t_pred_n];
end
if poly_degree >= 2
    A_poly_pred = [A_poly_pred, x_pred_n.^2, y_pred_n.^2, t_pred_n.^2, ...
                   x_pred_n.*y_pred_n, x_pred_n.*t_pred_n, y_pred_n.*t_pred_n];
end
if poly_degree >= 3
    A_poly_pred = [A_poly_pred, x_pred_n.^3, y_pred_n.^3, t_pred_n.^3, ...
                   x_pred_n.^2.*y_pred_n, x_pred_n.^2.*t_pred_n, ...
                   y_pred_n.^2.*x_pred_n, y_pred_n.^2.*t_pred_n, ...
                   t_pred_n.^2.*x_pred_n, t_pred_n.^2.*y_pred_n, ...
                   x_pred_n.*y_pred_n.*t_pred_n];
end

% Apply the same mask
A_poly_pred_red = A_poly_pred(:, mask_poly);

% Propagate to grid: diag(A * C_xx * A')
var_poly_grid = sum((A_poly_pred_red * C_xx_poly) .* A_poly_pred_red, 2);


% 5.2) Variance of the spline component
% 5.2.1) Rebuild the design matrix (A_obs) for the best model
% Needed to form the normal matrix (N) and get the covariance (C_cc)
n_obs_var = length(x_coords); 
A_obs_opt = sparse(n_obs_var, n_params);

% Pre-compute inverse deltas
inv_dX = 1/deltaX;
inv_dY = 1/deltaY;
inv_dT = 1/deltaT;

for n = 1:n_obs_var
    x = x_coords(n);
    y = y_coords(n);
    t = t_coords(n);
    
    % Basis indices
    i_x = floor((x - knots_x(1)) * inv_dX);
    if i_x >= num_x, i_x = num_x - 1; end
    if i_x < 0, i_x = 0; end
    
    i_y = floor((y - knots_y(1)) * inv_dY);
    if i_y >= num_y, i_y = num_y - 1; end
    if i_y < 0, i_y = 0; end
    
    i_t = floor((t - knots_t(1)) * inv_dT);
    if i_t >= num_t, i_t = num_t - 1; end
    if i_t < 0, i_t = 0; end
    
    % Local coords
    csi_x = (x - knots_x(i_x + 1)) * inv_dX;
    csi_y = (y - knots_y(i_y + 1)) * inv_dY;
    csi_t = (t - knots_t(i_t + 1)) * inv_dT;
    
    % Basis functions
    bx = bspline_basis(i_x, csi_x);
    by = bspline_basis(i_y, csi_y);
    bt = bspline_basis(i_t, csi_t);
    
    % Fill Row
    for l = -1:2
        idx_t = i_t + l + 1;
        for h = -1:2
            idx_y = i_y + h + 1;
            for k = -1:2
                idx_x = i_x + k + 1;
                if (idx_x >= 1 && idx_x <= n_ctrl_x && ...
                    idx_y >= 1 && idx_y <= n_ctrl_y && ...
                    idx_t >= 1 && idx_t <= n_ctrl_t)
                    
                    col = idx_y + (idx_x - 1) * n_ctrl_y + (idx_t - 1) * n_ctrl_y * n_ctrl_x;
                    A_obs_opt(n, col) = bx(k+2) * by(h+2) * bt(l+2);
                end
            end
        end
    end
end

% 5.2.2) Compute covariance matrix of coefficients (C_cc)
N_total = A_obs_opt' * A_obs_opt + lambda_spl * speye(n_params);

% Need the dense inverse for variance propagation
% C_cc = s02_spl * inv(N)
try
    C_cc = full(s02_spl * inv(N_total)); 
catch
    warning('Matrix close to singular. Using PINV for covariance.');
    C_cc = full(s02_spl * pinv(full(N_total)));
end
clear A_obs_opt N_total;

% 5.2.3) Propagate to prediction grid
x_vec_s = x_grid(:); % spatial points [n_spatial x 1]
y_vec_s = y_grid(:);
n_pixels = length(x_vec_s);
n_epochs_var = length(t_full);

var_spline_grid = zeros(n_spatial * n_epochs_var, 1);

for ti = 1:n_epochs_var
    curr_t = t_full(ti);
    
    % Time basis (constant for each epoch)
    i_t = floor((curr_t - knots_t(1)) * inv_dT);
    if i_t >= num_t, i_t = num_t - 1; end
    if i_t < 0, i_t = 0; end
    csi_t = (curr_t - knots_t(i_t + 1)) * inv_dT;
    bt = bspline_basis(i_t, csi_t);
    
    % Loop space
    for pi = 1:n_pixels
        curr_x = x_vec_s(pi);
        curr_y = y_vec_s(pi);
        
        % X basis
        i_x = floor((curr_x - knots_x(1)) * inv_dX);
        if i_x >= num_x, i_x = num_x - 1; end
        if i_x < 0, i_x = 0; end
        csi_x = (curr_x - knots_x(i_x + 1)) * inv_dX;
        bx = bspline_basis(i_x, csi_x);
        
        % Y basis
        i_y = floor((curr_y - knots_y(1)) * inv_dY);
        if i_y >= num_y, i_y = num_y - 1; end
        if i_y < 0, i_y = 0; end
        csi_y = (curr_y - knots_y(i_y + 1)) * inv_dY;
        by = bspline_basis(i_y, csi_y);
        
        % Collect indices and weights (vectorized dot product preparation)
        inds = zeros(64, 1);
        wts  = zeros(64, 1);
        count = 0;
        
        for l = -1:2
            idx_t = i_t + l + 1;
            for h = -1:2
                idx_y = i_y + h + 1;
                for k = -1:2
                    idx_x = i_x + k + 1;
                    if (idx_x >= 1 && idx_x <= n_ctrl_x && ...
                        idx_y >= 1 && idx_y <= n_ctrl_y && ...
                        idx_t >= 1 && idx_t <= n_ctrl_t)
                        
                        count = count + 1;
                        inds(count) = idx_y + (idx_x - 1) * n_ctrl_y + (idx_t - 1) * n_ctrl_y * n_ctrl_x;
                        wts(count)  = bx(k+2) * by(h+2) * bt(l+2);
                    end
                end
            end
        end
        
        % Calculate variance for this point: w' * C_cc * w
        if count > 0
            valid_inds = inds(1:count);
            valid_wts  = wts(1:count);
            
            % Extract sub-block of covariance
            cov_sub = C_cc(valid_inds, valid_inds);
            % Quadratic form
            var_val = valid_wts' * cov_sub * valid_wts;
            
            % Save to linear vector
            % Index corresponds to the logic: pixel varies fast, time varies slow
            global_idx = (ti - 1) * n_pixels + pi;
            var_spline_grid(global_idx) = var_val;
        end
    end
end


% 5.3) Total Variance and Reshape
% (assuming independence between poly and spline errors)
var_total_col = var_spline_grid;

% Reshape to [n_rows, n_cols, n_times]
final_std = reshape(sqrt(var_total_col), n_rows, n_cols, n_times);

disp('Variance propagation completed.');




%% 6) OUTPUT VARIABLES
disp('======== Step 6 ========');
disp('Output variables preparation started...');

% 6.1) Define latitude and longitude of output grid
[lat_full, lon_full] = utm2deg(x_grid(:), y_grid(:), repmat(utmZone, length(x_grid(:)), 1));
xyIN_AOI_flag = inpolygon(x_grid(:), y_grid(:), xyAOI(:,1), xyAOI(:,2));
lat_full_shp = lat_full(xyIN_AOI_flag);
lon_full_shp = lon_full(xyIN_AOI_flag);
lonlatIN_AOI_DET2D = [lon_full_shp, lat_full_shp];

% Lat/Lon coordinaates for PS plot
[lat_ps, lon_ps] = utm2deg(xyIN_AOI(:, 1), xyIN_AOI(:, 2), repmat(utmZone, size(xyIN_AOI, 1), 1));


% 6.2) Export final_signal_out for estimation grid points
% final_signal is [n_y x n_x x n_times], lonlatIN_AOI_STC2D is [n_x*n_y x 2]
% reshape final_signal to [n_spat x n_times], y-major order
n_spat = n_rows * n_cols;
final_signal_out = reshape(final_signal_shift, n_spat, n_times);  % [n_x*n_y x n_times]
final_std_out = reshape(final_std, n_spat, n_times);              % [n_x*n_y x n_times]

final_signal_out = final_signal_out(xyIN_AOI_flag, :);
final_std_out = final_std_out(xyIN_AOI_flag, :);

% verify sizes match
assert(size(final_signal_out, 1) == size(lonlatIN_AOI_DET2D, 1), ...
    'Mismatch: final_signal_out has %d rows, lonlatIN_AOI_STC2D has %d rows.', ...
    size(final_signal_out, 1), size(lonlatIN_AOI_DET2D, 1));
assert(size(final_signal_out, 2) == n_times, ...
    'Mismatch: final_signal_out has %d columns, expected %d times.', ...
    size(final_signal_out, 2), n_times);


% 6.3) Define the final time
dates_full = t_dateIN(1) + t_full;


% OUTPUT:
% lonlatIN_AOI_DET2D = lonlatIN_AOI_DET2D;
% dates_full = dates_full;
% t_full = t_full;
% final_signal_orig = final_signal_orig;
% final_signal = final_signal_out;
% final_std = final_std_out;

disp('Output variables preparation complete.');




%% 7) GIF WITH BASEMAP
disp('======== Step 7 ========');
disp('GIF creation started...');

% GIF creation
% R2026a-compat: getframe on invisible figure hangs, skip GIF on R2026a+
try
    gif_ok = isMATLABReleaseOlderThan('R2026a');  % getframe on invisible figure works
catch
    gif_ok = true;  % older MATLAB (pre-R2020b) lacks the check - assume OK
end
if gif_ok
figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
geoaxes;
v_min = min(round(prctile(final_signal_out(:), 5), 0), -5);
v_max = max(round(prctile(final_signal_out(:), 95), 0), 5);
h = waitbar(0, 'Plotting epochs...');
for t = 1:length(t_full)
    waitbar(t/length(t_full), h, sprintf('Plotting epoch %d/%d', t, length(t_full)));

    displ_at_t = final_signal_shift(:, :, t);
    displ_at_t = displ_at_t(:);
    displ_at_i_shp = displ_at_t(xyIN_AOI_flag);

    geobasemap satellite;
    if t == 1
        pause(5)
    end
    hold on;
    geoscatter(lat_full_shp, lon_full_shp, markerSize_DET2D, displ_at_i_shp, 'filled');
    geoscatter(lat_ps, lon_ps, markerSize_DET2D/4, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none');

    colormap(jet); clim([v_min, v_max]);
    c = colorbar; c.Label.String = 'LOS Displacement [mm]'; c.Label.FontSize = 15;
    title(sprintf('LOS Displacement on %s (2D)', datestr(dates_full(t))), 'FontSize', 18);
    hold off;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if t == 1
        imwrite(imind, cm, fullfile(figsDir, 'STdet_displ2D.gif'), 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, fullfile(figsDir, 'STdet_displ2D.gif'), 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
    clf;
end
close(h); close;
else
    disp('GIF creation skipped (R2026a+: getframe on invisible figure hangs)');
end % end if gif_ok


disp('GIF creation completed.');




%%             -------- END OF SCRIPT -------
disp('==================================================');
disp('--- 2D Deterministic-based modelling COMPLETED ---');
disp('==================================================');