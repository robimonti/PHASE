function result=nativeDet2D(state)
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
% Adapter guard: retain entries for skipped candidates before native sorting.
lambda_tmp=NaN(num_candidates,1);coef_tmp=cell(num_candidates,1);
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
    max_nnz = n_obs * 64; % max 64 non zero elements per row (4x4x4)
    I_idx = zeros(max_nnz, 1);
    J_idx = zeros(max_nnz, 1);
    V_val = zeros(max_nnz, 1);
    count_nnz = 0;
    
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

        % Fill Matrix in triplets
        for k = -1:2
            for h = -1:2
                for l = -1:2
                    abs_ix = i_x + k + 1; 
                    abs_iy = i_y + h + 1;
                    abs_it = i_t + l + 1;

                    if (abs_ix >= 1 && abs_ix <= n_ctrl_x && ...
                        abs_iy >= 1 && abs_iy <= n_ctrl_y && ...
                        abs_it >= 1 && abs_it <= n_ctrl_t)
                        
                        col = abs_iy + (abs_ix - 1) * n_ctrl_y + (abs_it - 1) * n_ctrl_y * n_ctrl_x;
                        
                        count_nnz = count_nnz + 1;
                        I_idx(count_nnz) = n;
                        J_idx(count_nnz) = col;
                        V_val(count_nnz) = alpha(k + 2, h + 2, l + 2);
                    end
                end
            end
        end
    end

    % Build the sparse matrix
    A = sparse(I_idx(1:count_nnz), J_idx(1:count_nnz), V_val(1:count_nnz), n_obs, n_params);

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
    N = A' * A + lambda * speye(n_params);
    rhs = A' * obs;
    if condest(N) > 1e10
        coef = pinv(full(N)) * rhs;
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



result=struct();
result.knots_x=knots_x;
result.knots_y=knots_y;
result.lambda_spl=lambda_spl;
result.s02_spl=s02_spl;
result.varNoise=varNoise;
result.knots_t=knots_t;
result.num_x=num_x;
result.num_y=num_y;
result.num_t=num_t;
result.coef_spl=coef_spl;
end
