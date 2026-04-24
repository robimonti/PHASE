function [lonlatIN_AOI_STC1D, dates_full, t_full, final_signal_orig, final_signal_out, final_std_out] = ...
    STmodel_STC1D(displIN_AOI, PSidIN_AOI, t_dateIN, t_relIN, centerline_data, figsDir, minMonths, ...
    gS_input_path, gS_output_path, gS_job_path, detectedOS, dtCov_STC1D, dsCov_STC1D, utmZone, xyIN_AOI, step_t, ...
    detrend_method, use_inclined_means, poly_degree, tCovModel_STC1D, sCovModel_STC1D, markerSize_STC1D)

% STmodel_STC1D Performs spatio-temporal stochastic modelling for displacement data
%
% Inputs:
%   displIN_AOI           - n x m matrix of displacement data (n points, m epochs)
%   PSidIN_AOI            - n x 1 vector of point IDs
%   t_dateIN              - m x 1 vector of absolute dates (datetime)
%   t_relIN               - m x 1 vector of relative times (days)
%   centerline_data       - Structure with fields:
%                           - xy_centerline: k x 2 matrix of centerline coordinates
%                           - s: k x 1 vector of cumulative distances along centerline
%                           - projected_distances: n x 1 vector of projected distances
%   figsDir               - Directory path for saving figures
%   file_out              - Minimum period length for splines interpolation
%   gS_input_path         - path to input folder for geoSplinter
%   gS_output_path        - path to output folder for geoSplinter
%   gS_job_path           - path to job folder for geoSplinter
%   detectedOS            - detected operative system
%   dtCov_STC1D           - step size to average time covariance (NaN means automatic)
%   dsCov_STC1D           - step size to average space covariance (NaN means automatic)
%   utmZone               - UTM zone for reprojection of coordinates
%   xyIN_AOI              - coordinates of PS inside the AOI
%   step_t                - time step for the estimation
%   detrend_method        - option to select the processing method: residual atmospheric errors (residualAtm) or clean data (cleanObs)
%   use_inclined_means    - flag for using flat or inclines plane for residual atmospheric removal
%   poly_degree           - maximum polynomial degree for the cleanObs case
%   tCovModel_STC2D       - temporal covariance model
%   sCovModel_STC2D       - spatial covariance model
%   markerSize_STC2D      - marker size for the GIF plot
%
% Output:
%   lonlatIN_AOI_STC1D    - lon/lat coordinates of estimated centerline points
%   dates_full            - dates prediction
%   t_full                - relative time prediction
%   final_signal_orig     - modelled signal in observation coordinates
%   final_signal_out      - modelled signal in query coordinates
%   final_std_out         - uncertainty
%
%
% !! WARNING: This function can be computationally intensive on your RAM,
%             so verify it. In case of problems try reducing the resolution
%             of your output from the configuration GUI.


% validate inputs
if nargin < 22
    error('STC 1D: All input arguments are required.');
end
% validate centerline inputs
if isempty(centerline_data) || isempty(centerline_data.xy_centerline) || isempty(centerline_data.s) || isempty(centerline_data.projected_distances)
    error('STC 1D: centerline_data or its fields are empty for 1D interpolation.');
end


disp('--- 1D Stochastic-based modelling started... ---');


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

% Create the full grid of observed positions (x, y) and times
grid_t0 = repmat(t_relIN, np, 1);  % [np x nt]
grid_s0 = repmat(s_SAR, 1, nt);  % [np x nt]

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




%% 3) LEAST SQUARES COLLOCATION
disp('======== Step 3 ========');
disp('Iterative 2D Least Squares Collocation started...');

% observations (start with polynomial residuals)
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

iter = 0;
delta_conv = 1000;
spatialMeans_lsc_itPrec = Inf(length(t_relIN), 1);

% For parallel computing
maxLen = max(sum(~isnan(dataRawCleanT_sorted), 2));
nWorkers = autoNumWorkers(maxLen);

% Final outputs (overwritten only at the last iteration)
full_poly_grid      = zeros(size(residuals_grid_full));
spatialMeans_grid   = zeros(size(residuals_grid_full));
C_obs_obs_m         = [];

while iter < 2
    iter = iter + 1;
    fprintf('\n=== Full Iterative LSC - Iteration %d ===\n', iter);

    % 3.1) Covariance matrix Q_ll (same for trend and spatial means)
    if iter == 1
        variances = ones(length(y0_lsc), 1);
        idx_master = (grid_t_flt_col == master_t_rel);
        if any(idx_master)
            variances(idx_master) = 0.1;
        end
        Q_ll = spdiags(variances, 0, length(y0_lsc), length(y0_lsc));
    else
        Q_ll = C_obs_obs_m;  % <- already contains both signal and noise 
    end


    % 3.2) Design matrix construction
    n_epochs = length(t_relIN);

    if strcmp(detrend_method, 'residualAtm')

        % 3.2.A) Residual errors in the data (Orbits + Atmosphere)

        % Global orbital trend (static plane a*x + b*y)
        A_trend = s_centered;
        n_trend_cols = size(A_trend, 2);
        
        if ~use_inclined_means
            % 3.2.A1 Small area: constant atmosphere per epoch
           
            n_PS = size(residuals_grid_full, 1);
            
            % Identity matrix for epochs, replicated for each PS
            A_atmos_full = repmat(eye(n_epochs), n_PS, 1);
            
            % Extract valid observations
            A_atmos = A_atmos_full(~mask_obs_col, :);
            
        else
            % 3.2.A2) Large area: ramped atmosphere per epoch
            
            % 2 parameters per epoch: [Bias, SlopeS]
            % Total atmos params = 2 * n_epochs
            
            % For each observation (i, t) -> columns corresponding to epoch t.
            % The row is: [ ... 0 0 1 s_i 0 0 ... ] at the columns for epoch t.
            
            % Map each valid observation to its epoch index (1 to n_epochs)
            [~, obs_epoch_idx] = ismember(grid_t_flt_col, t_relIN);
            
            % Base indices for the 3 parameters of each epoch
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

    elseif strcmp(detrend_method, 'cleanObs')

        % 3.2.B) Polynomial surface (empirical)

        % Normalize coordinates for numerical stability
        % (using the vectors corresponding to current VALID observations)
        t_n = t_centered;
        s_n = s_centered;
        
        % Start with intercept
        A_poly = ones(length(y0_lsc), 1);
        
        % Degree 1
        if poly_degree >= 1
            A_poly = [A_poly, s_n, t_n];
        end
        
        % Degree 2
        if poly_degree >= 2
            A_poly = [A_poly, s_n.^2, t_n.^2, s_n.*t_n];
        end
        
        % Degree 3
        if poly_degree >= 3
            A_poly = [A_poly, s_n.^3, t_n.^3, ...
                      s_n.^2.*t_n, t_n.^2.*s_n];
        end
        
        % ->
        A_full = A_poly;
        n_trend_cols = size(A_full, 2);

    end


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

    if strcmp(detrend_method, 'residualAtm')
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
            spatialMeans_grid_temp = repmat(x_atmos_all', size(residuals_grid_full, 1), 1);
            
        else
            % Case: Inclined Planes (Ramps)
            % x_atmos_all is [2*n_epochs x 1]
            % Reshape to [2 x n_epochs] -> Row 1: Bias, Row 2: SlopeS
            atmos_params = reshape(x_atmos_all, 2, n_epochs);
            
            biases  = atmos_params(1, :); % [1 x n_epochs]
            slopesS = atmos_params(2, :);
            
            % Construct the full grid correction
            % Atmos(t,x,y) = Bias(t) + Ss(t)*s
            
            % Grid X and Y (centered)
            S_grid = reshape(s_full_centered, length(t_relIN), [])'; % [n_PS x n_epochs]
            
            % Expand parameters to grid size
            B_grid = repmat(biases, size(residuals_grid_full, 1), 1);
            Ss_grid = repmat(slopesS, size(residuals_grid_full, 1), 1);
            
            spatialMeans_grid_temp = B_grid + (Ss_grid .* S_grid);
        end

        x_means = x_atmos_all;
        
    elseif strcmp(detrend_method, 'cleanObs')
        x_poly_coeffs = x_est_all;      % all parameters are poly coeffs
        x_means = zeros(n_epochs, 1);   % no separate means
        
        % Build full grid A_poly
        t_full_n = t_full_centered;
        s_full_n = s_full_centered;
        
        A_full_grid = ones(length(t_full_n), 1);

        if poly_degree >= 1
            A_full_grid = [A_full_grid, s_full_n, t_full_n];
        end
        
        if poly_degree >= 2
             A_full_grid = [A_full_grid, s_full_n.^2, t_full_n.^2, s_full_n.*t_full_n];
        end
        
        if poly_degree >= 3
             A_full_grid = [A_full_grid, s_full_n.^3, t_full_n.^3, ...
                      s_full_n.^2.*t_full_n, t_full_n.^2.*s_full_n];
        end
        
        trend_full_col = A_full_grid * x_poly_coeffs;
        trend_grid_full = reshape(trend_full_col, length(t_relIN), [])';
        
        spatialMeans_grid_temp = zeros(size(trend_grid_full));   % zero
    end
    
    fprintf('    > Elimination Complete. %d parameters remaining.\n', sum(mask_params));


    % 3.5) Remove both components from current residuals
    residuals_grid_full = dataRawCleanT_sorted ...
        - trend_grid_full ...
        - spatialMeans_grid_temp;


    % 3.6) Convergence check
    delta_conv = max(abs(x_means - spatialMeans_lsc_itPrec));
    fprintf('Max change in spatial means: %.4f mm\n', delta_conv);
    spatialMeans_lsc_itPrec = x_means;


    % 3.7) Global Empirical Covariance Estimation
    disp('Computing global empirical covariance...');
    
        % 3.7.1) Covariance in TIME
        if isnan(dtCov_STC1D)
            dtCov = mode(diff(t_relIN)); if dtCov==0, dtCov=1; end
        else
            dtCov = dtCov_STC1D;
        end
        max_lag_time = ceil((t_relIN(end)-t_relIN(1))/dtCov) + 1;
        
        num_ps = size(residuals_grid_full, 1);
        accum_T_cells = cell(num_ps, 1);
        var_T_cells   = cell(num_ps, 1);   % store variance parts
        
        parfor i = 1:num_ps
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
        Variance_T = Var_Global_T(1) / Var_Global_T(2);
        
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
        if isnan(dsCov_STC1D)
            dsCov = mode(diff(s_SAR)); if dsCov==0, dsCov=1; end
        else
            dsCov = dsCov_STC1D;
        end 
        edges = 0:dsCov:ceil(max(s_SAR) / dsCov) * dsCov;
        n_epochs = size(residuals_grid_full, 2);
        accum_S_cells = cell(n_epochs, 1);
        var_S_cells   = cell(n_epochs, 1);
        var_S_cells1   = cell(n_epochs, 1);
    
        % Broadcast variable for parfor
        method_flag = detrend_method; 
        
        parfor i = 1:n_epochs
            sp = residuals_grid_full(:, i);
            valid = ~isnan(sp);
            
            if sum(valid) > 1
                sp_val = sp(valid);
                coords_val = s_SAR(valid);
                
                % Conditional centering
                switch method_flag
                    case 'cleanObs'
                    % In cleanObs mode, the global surface removal might leave 
                    % epoch-wise offsets. The mean is removed locally strictly 
                    % for covariance estimation to ensure the tail decays to zero.
                    sp_input = sp_val - mean(sp_val);
    
                    [~, s_sq1, n_p1] = f1DEmpCovEst_accum(sp_val, coords_val, edges, true);
                    var_S_cells1{i}   = [s_sq1, n_p1];
                    case 'residualAtm'
                    % In residualAtm mode, the iterative LSC has already removed 
                    % the atmosphere (mean/ramp).
                    sp_input = sp_val;
                end
                
                [local_acc, s_sq, n_p] = f1DEmpCovEst_accum(sp_input, coords_val, edges, true);
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
        Variance_S = Var_Global_S(1) / Var_Global_S(2);
        
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
        lag_t1 = 0:round(lag_t(end)/100):round(lag_t(end));
    else
        lag_t1 = lag_t;
    end
    y_data_t = eCovF_mean_t;
    w_data_t = Cecf_mean_t;
    
    % Use the chosen model
    switch tCovModel_STC1D
        case 'exponential'
            [c_opt_t, mCovF1] = fit_exp_auto(x_data_t, y_data_t, w_data_t);
        case 'gaussian'
            [c_opt_t, mCovF1] = fit_gaussian_auto(x_data_t, y_data_t, w_data_t);
        case 'gaussian_cos'
            [c_opt_t, mCovF1] = fit_gaussian_cosine_auto(x_data_t, y_data_t, w_data_t);
    end
    c1_t = c_opt_t;
    
    
    % 3.9) Fit SPACE covariance
    x_data_s = lag_s;
    if length(lag_s) < 20
        lag_s1 = 0:round(lag_s(end)/100):round(lag_s(end));
    else
        lag_s1 = lag_s;
    end
    y_data_s = eCovF_mean_s;
    w_data_s = Cecf_mean_s;
    
    % Use the chosen model
    switch sCovModel_STC1D
        case 'exponential'
            [c_opt_s, mCovF2] = fit_exp_auto(x_data_s, y_data_s, w_data_s);
        case 'gaussian'
            [c_opt_s, mCovF2] = fit_gaussian_auto(x_data_s, y_data_s, w_data_s);
        case 'gaussian_cos'
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

    % Compute the covariance matrix of observations
    xy_obs_all = [grid_t_flt_col, grid_s_flt_col];
    [C_obs_obs_m, ~, C_vv_m] = compute_C_obs_obs_nugget(xy_obs_all, ...
        mCovF1, c1_t, mCovF2, c1_s, A_noise, master_t_rel);
        % → This updates the covariance for the next iteration

    % 3.11) Save final results only on the last iteration
    if delta_conv < 1
        full_poly_grid     = trend_grid_full;
        spatialMeans_grid  = spatialMeans_grid_temp;
        % C_obs_obs_m already contains the final one from above
        fprintf('Converged after %d iterations!\n', iter);
    end

end


% Final cleaned time series for collocation
data_final_clean = dataRawCleanT_sorted ...
    - full_poly_grid ...
    - spatialMeans_grid;

delete(gcp('nocreate'));


% --- Plot covariances ---
figure('Visible', 'off', 'Position', [100, 100, 1400, 600]);

ax1 = subplot(1,2,1);
scatter(lag_t, eCovF_mean_t, 30, 'k', 'filled')
hold on; plot(lag_t1, mCovF1(c1_t, lag_t1), 'm', 'LineWidth', 2)
plot([1, 1], [mCovF1(c1_t, 0), eCovF_mean_t(1)], 'r--', 'LineWidth', 2)
xlim([0 round(lag_t(end)/2)]); grid on;
legend('empirical', 'signal', 'noise')
xlabel('Lags [days]'); ylabel('Covariance [mm^2]')
set(gca, 'FontSize', 15)

ax2 = subplot(1,2,2);
scatter(lag_s, eCovF_mean_s, 30, 'k', 'filled', 'DisplayName', 'empirical')
hold on; plot(lag_s1, mCovF2(c1_s, lag_s1), 'm', 'LineWidth', 2, 'DisplayName', 'signal')
plot([0.5, 0.5], [mCovF2(c1_s, 0), eCovF_mean_s(1)], 'r--', 'LineWidth', 2, 'DisplayName', 'noise');
xlim([0 round(lag_s(end)/2)]); grid on;
xlabel('Lags [m]'); ylabel('Covariance [mm^2]'); legend show
set(gca, 'FontSize', 15)

yl = [ ...
    min([ylim(ax1), ylim(ax2)]), ...
    max([ylim(ax1), ylim(ax2)]) ];
ylim(ax1, yl)
ylim(ax2, yl)

print(gcf, fullfile(figsDir, 'STC_cov2D.png'), '-dpng', '-r300')


% --- end ---

disp('Iterative 2D Least Squares Collocation completed.');




%% 4) COLLOCATION - PREDICTION
disp('======== Step 4 ========');
disp('Collocation prediction started...');

% 4.1) Prepare observations (xy_obs)
xy_obs = [grid_t_flt_col, grid_s_flt_col];
valid_obs = ~isnan(xy_obs(:, 1)) & ~isnan(xy_obs(:, 2));
xy_obs = xy_obs(valid_obs, :);

% Prepare data vector (residuals from LSC loop)
data_final_clean_col = reshape(data_final_clean', [], 1);
data_final_clean_col(isnan(data_final_clean_col)) = [];
data_final_clean_col = data_final_clean_col(valid_obs);

if length(data_final_clean_col) ~= sum(valid_obs)
    error('Mismatch: data length (%d) does not match valid_obs (%d)', ...
          length(data_final_clean_col), sum(valid_obs));
end


% 4.2) Prepare estimation grid (xy_est)
% create the full grid of estimation distances and times
t_full = (t_relIN(1):step_t:t_relIN(end))'; % [n_t x 1]

% Define spatial grid 
s_full = s_full_SAR; 

grid_t_full = repmat(t_full', length(s_full), 1);
grid_s_full = repmat(s_full, 1, length(t_full)); 

% Flatten to column vectors (time-major order)
grid_t_full_col = reshape(grid_t_full, [], 1); 
grid_s_full_col = reshape(grid_s_full, [], 1);

xy_est = [grid_t_full_col, grid_s_full_col]; 
clear grid_t_full;


% 4.3) Compute cross-covariance matrices (C_est_obs)
% Observations covariance matrix (from LSC loop)
C_obs_obs = C_obs_obs_m;
clear C_obs_obs_m

% Extract coordinates
t_est = xy_est(:, 1);
s_est_vec = xy_est(:, 2);

t_obs = xy_obs(:, 1);
s_obs = xy_obs(:, 2);

% 4.3.1) Temporal distance (tau_t)
tau_t = abs(t_est - t_obs'); 

% 4.3.2) Spatial distance (tau_s)
tau_s = abs(s_est_vec - s_obs');

% 4.3.3) Evaluate covariance functions
Ct = mCovF1(c1_t, tau_t);         
Cs_signal = mCovF2(c1_s, tau_s);  

% 4.3.4) Combine (separable model)
C_est_obs = (Ct .* Cs_signal) ./ c1_t(1);


% 4.4) Solve Collocation (stochastic component)
weights = C_obs_obs \ data_final_clean_col;
coll_est = C_est_obs * weights; 

% Reshape stochastic result to [n_s x n_t]
n_s = length(s_full);
n_t = length(t_full);

coll_est_2D = reshape(coll_est, n_s, n_t); 


% 4.5) Reconstruct deterministic trend on estimation grid - only for cleanObs case -
% Centering (must match the loop logic)
s_est_centered = s_est_vec - s0;
t_est_centered = t_est - t0; 

switch detrend_method
    case 'cleanObs'
    % model: poly(t,x)
    % params: x_poly_coeffs
    
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

    % reshape trend to grid
    poly_fit_full_2D = reshape(z_poly_est, n_s, n_t);

end


% 4.6) Final summation
switch detrend_method
    case 'residualAtm'
    % The trend (x_trend) was defined as "Orbital Error".
    % The means (x_means) were defined as "Atmospheric Artifacts".
    % -> The final signal is just the stochastic deformation.
    
    final_signal = coll_est_2D;

    poly_fit_full_2D = zeros(size(final_signal));

    case 'cleanObs'
    % The trend (x_poly_coeffs) was defined as the "deterministic signal".
    % Restore it to the stochastic deformation to get the total field.
    
    final_signal = coll_est_2D + poly_fit_full_2D;

end


% 4.7) Apply reference shift (relative to first epoch, t=0)
% This ensures the deformation starts at 0 for all pixels.
final_signal_shift = final_signal - final_signal(:,1);

% Also shift the components for plotting consistency
coll_est_grid_shift = coll_est_2D - coll_est_2D(:,1);
poly_fit_full_grid_shift = poly_fit_full_2D - poly_fit_full_2D(:,1);


% 4.8) Prediction on observation points
% Define target coordinates (full original grid)
% These variables were created in Step 2
target_t = t_full_0; % [np*nt x 1]
target_s = s_full_0;

% Source coordinates (valid data used for model)
% t_obs, s_obs are already defined in Step 4.1

% 4.8.1) Compute stochastic component on full grid
% Temporal distance
tau_t_full = abs(target_t - t_obs'); % [np*nt x n_valid_obs]

% Spatial distance
tau_s_full = abs(target_s - s_obs');

% Evaluate signal covariance (no noise)
% C_target_source(i,j) = Ct(tau) * Cs(tau)
C_target_source = (mCovF1(c1_t, tau_t_full) .* mCovF2(c1_s, tau_s_full)) ./ c1_t(1);

% Stochastic prediction: s = C_target_source * weights
% weights = inv(C_obs_obs) * residuals (calculated in 4.4)
signal_stoch_full = C_target_source * weights;

% 4.8.2) Compute deterministic component at observation points
signal_det_full = zeros(size(signal_stoch_full));

switch detrend_method
    case 'cleanObs'
        % Reconstruct the polynomial trend for observations
        % centering (as Step 3)
        t_target_n = target_t - t0;
        s_target_n = target_s - s0;
    
        % Build design matrix (A_poly_target)
        A_poly_target = ones(length(t_target_n), 1);

        if poly_degree >= 1
            A_poly_target = [A_poly_target, s_target_n, t_target_n];
        end
    
        if poly_degree == 2
            A_poly_target = [A_poly_target, s_target_n.^2, t_target_n.^2, ...
                          s_target_n.*t_target_n];
        end
        if poly_degree == 3
            A_poly_target = [A_poly_target, s_target_n.^3, t_target_n.^3, ...
                          s_target_n.^2.*t_target_n, t_target_n.^2.*s_target_n];
        end
    
        % Apply coefficients (x_poly_coeffs from Step 3)
        signal_det_full = A_poly_target * x_poly_coeffs;

    case 'residualAtm'
        % The desired signal is purely the stochastic deformation.
        signal_det_full = zeros(size(signal_stoch_full));

end

% 4.8.3) Final summation and reshape
final_signal_orig_col = signal_stoch_full + signal_det_full;
final_signal_orig = reshape(final_signal_orig_col, nt, np)';

disp('Collocation prediction completed.');




%% 5) TOTAL VARIANCE PROPAGATION
disp('======== Step 5 ========');
disp('Variance propagation started...');

% 5.1) Setup target coordinates
% xy_est is [n_total x 2] -> [Time, S]
n_2D_total = size(xy_est, 1); 


% 5.2) Construct design matrices (A_est_2D and A_obs)
switch detrend_method
    case 'residualAtm'
    % - Global orbital trend
    % Prediction grid
    A_2D_trend = xy_est(:,2) - s0;
    % Observation grid
    A_obs_trend = xy_obs(:,2) - s0;
    
    % - Atmosphere (means or ramps)
    unique_epochs = t_relIN; 
    n_epochs_LSC = length(unique_epochs);
    n_obs = size(y0_lsc, 1);
    
    if ~use_inclined_means
        % Small AOI: constant means (1 param per epoch)
        n_atmos_params = n_epochs_LSC;
        
        % Prediction: atmosphere is NOT predicted (zeros)
        A_2D_atmos = zeros(n_2D_total, n_atmos_params);
        
        % Observations: identity mapping
        A_obs_atmos = zeros(n_obs, n_atmos_params);
        [~, epoch_idx] = ismember(xy_obs(:,1), unique_epochs);
        valid_epochs = epoch_idx > 0;
        linear_idx = sub2ind(size(A_obs_atmos), find(valid_epochs), epoch_idx(valid_epochs));
        A_obs_atmos(linear_idx) = 1;
        
    else
        % Large AOI: inclined means (2 params per epoch)
        n_atmos_params = 2 * n_epochs_LSC;
        
        % Prediction: atmosphere is NOT predicted (zeros)
        A_2D_atmos = zeros(n_2D_total, n_atmos_params);
        
        % Observations: block ramps (Bias, Ss)
        [~, obs_epoch_idx] = ismember(xy_obs(:,1), unique_epochs);
        
        % Column indices
        col_bias   = (obs_epoch_idx - 1) * 3 + 1;
        col_slopeS = (obs_epoch_idx - 1) * 2 + 2;
        
        n_curr_obs = size(xy_obs, 1);
        
        % Build sparse matrix
        I = repmat((1:n_curr_obs)', 2, 1);
        J = [col_bias; col_slopeS];
        % Values: 1, (x-x0)
        V = [ones(n_curr_obs, 1); (xy_obs(:,2)-s0)];
        
        A_obs_atmos = sparse(I, J, V, n_curr_obs, n_atmos_params);
        if n_curr_obs < 10000, A_obs_atmos = full(A_obs_atmos); end
    end
    
    % Combine
    A_est_2D = [A_2D_trend, A_2D_atmos];
    A_obs     = [A_obs_trend, A_obs_atmos];

    case 'cleanObs'
    % Prepare prediction coordinates
    t_2D_n = xy_est(:,1) - t0;
    s_2D_n = xy_est(:,2) - s0;
    
    % Prepare observation coordinates
    t_obs_n = xy_obs(:,1) - t0;
    s_obs_n = xy_obs(:,2) - s0;
    
    % - Build A_poly_2D
    A_poly_2D = ones(length(t_2D_n), 1);

    if poly_degree >= 1
        A_poly_2D = [A_poly_2D, s_2D_n, t_2D_n];
    end
    
    if poly_degree >= 2
        A_poly_2D = [A_poly_2D, s_2D_n.^2, t_2D_n.^2, s_2D_n.*t_2D_n];
    end
    if poly_degree >= 3
        A_poly_2D = [A_poly_2D, s_2D_n.^3, t_2D_n.^3, ...
                      s_2D_n.^2.*t_2D_n, t_2D_n.^2.*s_2D_n];
    end
    
    % - Build A_poly_obs
    A_poly_obs = ones(length(t_obs_n), 1);

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
    
    A_est_2D = A_poly_2D;
    A_obs     = A_poly_obs;
end

% Check compatibility with mask
if size(A_est_2D, 2) ~= length(mask_params)
    error('Design matrix width (%d) does not match mask_params (%d).', ...
        size(A_est_2D, 2), length(mask_params));
end

% Filter based on backward elimination results
A_est_2D = A_est_2D(:, mask_params);
A_obs     = A_obs(:, mask_params);


% 5.3) Construct C_est_obs_2D
% Temporal distance
tau_t_2D = abs(xy_est(:,1) - xy_obs(:,1)');

% Spatial Euclidean distance
tau_s_2D = abs(xy_est(:,2) - xy_obs(:,2)');

% Combine
C_est_obs_2D = (mCovF1(c1_t, tau_t_2D) .* mCovF2(c1_s, tau_s_2D)) ./ c1_t(1);

% 5.4) Signal covariance (diagonal of C_sig)
C_sig_diag = repmat(c1_s(1), n_2D_total, 1);


% 5.5) Calculation with Cholesky decomposition method
% Compute L
L = chol(C_obs_obs, 'lower'); 

% 5.5.1) Collocation weight term: Y = L \ C_est_obs_2D'
Y = L \ C_est_obs_2D'; 
term_coll_reduction = sum(Y.^2, 1)'; % [n_2D x 1]

% 5.5.2) Trend weight term: Z = L \ A_obs
Z = L \ A_obs;
N_mat = Z' * Z; 
N_inv = inv(N_mat);

% 5.5.3) Propagation factor: F = A_est_2D - (Y' * Z)
F = A_est_2D - (Y' * Z);

% 5.5.4) Trend variance
term_trend_variance = sum((F * N_inv) .* F, 2);


% 5.6) Final recombination
var_2D = C_sig_diag - term_coll_reduction + term_trend_variance;

    % Sanity check
    var_2D(var_2D < 0) = 0;

% Reshape for output
final_std = reshape(sqrt(var_2D), n_s, n_t);
final_std_noise_3D = reshape(sqrt(var_2D + A_noise), n_s, n_t);

disp('Variance propagation completed.');




%% 6) OUTPUT VARIABLES
disp('======== Step 6 ========');
disp('Output variables preparation started...');

% 6.1) Define latitude and longitude of output grid
grid_x = centerline_data.xy_centerline(:,1);
grid_y = centerline_data.xy_centerline(:,2);
[lat_centerline, lon_centerline] = utm2deg(grid_x, grid_y, repmat(utmZone, length(grid_x), 1));
lonlatIN_AOI_STC1D = [lon_centerline, lat_centerline];

% Lat/Lon coordinaates for PS plot
[lat_ps, lon_ps] = utm2deg(xyIN_AOI(:, 1), xyIN_AOI(:, 2), repmat(utmZone, size(xyIN_AOI, 1), 1));


% 6.2) Export final_signal_out for estimation grid points
final_signal_out = final_signal_shift;
final_std_out = final_std;

% verify sizes match
assert(size(final_signal_out, 1) == size(lonlatIN_AOI_STC1D, 1), ...
    'Mismatch: final_signal_out has %d rows, lonlatIN_AOI_STC2D has %d rows.', ...
    size(final_signal_out, 1), size(lonlatIN_AOI_STC1D, 1));
assert(size(final_signal_out, 2) == n_t, ...
    'Mismatch: final_signal_out has %d columns, expected %d times.', ...
    size(final_signal_out, 2), n_t);


% 6.3) Define the final time
dates_full = t_dateIN(1) + t_full;


% 6.4) Re-sort export variable final_signal_orig to match the input one
[~, unsort_idx] = sort(idx_SAR);
final_signal_orig = final_signal_orig(unsort_idx, :);


% OUTPUT:
% lonlatIN_AOI_STC2D = lonlatIN_AOI_STC1D;
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

    displ_at_t = final_signal_shift(:, t);
    displ_at_t_shp = displ_at_t(:);

    geobasemap satellite;
    if t == 1
        pause(20)
    end
    hold on;
    geoscatter(lat_centerline, lon_centerline, markerSize_STC1D, displ_at_t_shp, 'filled');
    geoscatter(lat_ps, lon_ps, markerSize_STC1D/4, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none');

    colormap(jet); clim([v_min, v_max]);
    c = colorbar; c.Label.String = 'LOS Displacement [mm]'; c.Label.FontSize = 15;
    title(sprintf('LOS Displacement on %s (1D)', datestr(dates_full(t))), 'FontSize', 18);
    hold off;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if t == 1
        imwrite(imind, cm, fullfile(figsDir, 'STstc_displ1D.gif'), 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, fullfile(figsDir, 'STstc_displ1D.gif'), 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
    clf;
end
close(h); close;
else
    disp('GIF creation skipped (R2026a+: getframe on invisible figure hangs)');
end % end if gif_ok



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
filename = strcat(figsDir, filesep, 'STstc_displ1D.png');
print(gcf, filename, '-dpng', '-r300')
close

disp('GIF creation completed.');




%%           -------- END OF SCRIPT -------
disp('===============================================');
disp('--- 1D Stochastic-based modelling COMPLETED ---');
disp('===============================================');