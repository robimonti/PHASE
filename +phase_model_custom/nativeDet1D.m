function result=nativeDet1D(state)
% Original residual spline selection/fit, copied verbatim; no model reduction.
names=fieldnames(state);
for variableIndex=1:numel(names)
    eval([names{variableIndex} '=state.(names{variableIndex});']);
end
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
        gS_exec = fullfile(gS_dir, 'geoSplinter_analysis');
        gS_job_file = fullfile(gS_job_path, [file_out, '.job']);
        [status, geoOutput] = runGeoSplinter(gS_exec,gS_job_file);
        if status ~= 0
            error('Error executing geoSplinter_analysis for file %s: %s', ...
                file_out,strtrim(geoOutput));
        end

        % Import results
            %[~, ~, ~,~] = geoSplinter(fullfile(gS_output_path, file_out), ext);
        [gS_out_data_tmp, gS_out_raster, ~, ~, gS_out_nodes] = geoSplinter_noFig(...
            fullfile(gS_output_path, file_out), ext);

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
        n_points_spl = length(x_coords_sort);
        max_nnz = n_points_spl * 16; % 4x4 spline support
        I_idx = zeros(max_nnz, 1);
        J_idx = zeros(max_nnz, 1);
        V_val = zeros(max_nnz, 1);
        count_nnz = 0;
        
        % Compute design matrix
        for n = 1:n_points_spl
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
                            
                            col = (i_y + h) + (i_x + k) * length(knots_y) + 1;
                            
                            count_nnz = count_nnz + 1;
                            I_idx(count_nnz) = n;
                            J_idx(count_nnz) = col;
                            V_val(count_nnz) = alpha(k + 2, h + 2);
                        end
                    end
                end
            end
        end

        % Create the spare matrix
        A_spline = sparse(I_idx(1:count_nnz), J_idx(1:count_nnz), V_val(1:count_nnz), n_points_spl, n_params);

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



result=struct();
result.knots_x=knots_x;
result.knots_y=knots_y;
result.lambda_spl=lambda_spl;
result.s02_spl=s02_spl;
result.varNoise=varNoise;
result.row_col_pairs=row_col_pairs;
result.var_idx=var_idx;
result.gS_output_path=gS_output_path;
pair=row_col_pairs(var_idx,:); file=sprintf('DET1D_bic_r%d_c%d',pair(1),pair(2));
result.nativeCoefficients=load(fullfile(gS_output_path,[file '.par.txt']));
end
