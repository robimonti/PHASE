function [C_obs_obs_fin, C_obs_obs, Cvv_obs] = compute_C_obs_obs_nugget_2D(xy_obs, mCovF1, c1_t, mCovF2, c1_s, A_noise, t_master)
 
% compute_C_obs_obs_nugget_2D Computes the combined covariance matrix C_obs_obs for 2D Spatial + Time data.
%
% Inputs:
%   xy_obs   - n_obs x 3 matrix: [Time, X, Y]
%   mCovF1   - Function handle for temporal covariance
%   c1_t     - Parameters for temporal covariance
%   mCovF2   - Function handle for spatial covariance
%   c1_s     - Parameters for spatial covariance
%   A_noise  - Amplitude for spatial noise (Nugget)
%   t_master - Master epoch time (for variance reduction)
%
% Output:
%   C_obs_obs_fin - Signal + Noise matrix
%   C_obs_obs     - Signal only matrix
%   Cvv_obs       - Noise only matrix

% Number of observations
n_obs = size(xy_obs, 1);

if size(xy_obs, 2) < 3
    error('xy_obs must have 3 columns: [Time, X, Y]');
end

% 1. Extract Vectors
t_vec = xy_obs(:, 1);
x_vec = xy_obs(:, 2);
y_vec = xy_obs(:, 3);

% 2. Calculate Distance Matrices (Vectorized)

% Temporal Distance (Tau_t = |t_i - t_j|)
Tau_t = abs(t_vec - t_vec');

% Spatial Euclidean Distance (Tau_s = sqrt((x_i-x_j)^2 + (y_i-y_j)^2))
Diff_X = x_vec - x_vec';
Diff_Y = y_vec - y_vec';
Tau_s  = sqrt(Diff_X.^2 + Diff_Y.^2);

% 3. Evaluate Covariance Functions
% Temporal Component
Ct_mat = mCovF1(c1_t, Tau_t);

% Spatial Component
Cs_mat = mCovF2(c1_s, Tau_s);

% 4. Combine (Separable Model)
% C(t,s) = C(t) * C(s) / Variance_Time 
% (Standard separable assumption: C_space * C_time normalized)
if c1_t(1) ~= 0
    C_obs_obs = (Ct_mat .* Cs_mat) ./ c1_t(1);
else
    warning('Temporal variance is zero. Returning zero covariance.');
    C_obs_obs = zeros(n_obs);
end

% 5. Construct Noise Matrix (Diagonal)
% Create a vector of noise variances
diag_noise = repmat(A_noise, n_obs, 1);

% Reduce noise for master epoch (constraint enforcement)
% Using a small tolerance for float comparison
idx_master = abs(t_vec - t_master) < 1e-5;
if any(idx_master)
    diag_noise(idx_master) = A_noise * 0.1;
end

% Create diagonal matrix
Cvv_obs = diag(diag_noise);

% 6. Final Sum
C_obs_obs_fin = C_obs_obs + Cvv_obs;

% 7. Eigenvalue check (Optional Debugging)
% if n_obs <= 2000
%     e = eig(C_obs_obs_fin);
%     if min(e) < 0
%         fprintf('Warning: Covariance matrix is not positive definite.\n');
%         fprintf('Min eigenvalue   : %g\n', min(e));
%     end
% end
end