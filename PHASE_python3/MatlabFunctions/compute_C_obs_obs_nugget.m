function [C_obs_obs_fin, C_obs_obs, Cvv_obs] = compute_C_obs_obs_nugget(xy_obs, mCovF1, c1_t, mCovF2, c1_s, A_noise, t_master)

% compute_C_obs_obs_nugget Computes the combined covariance matrix C_obs_obs
%
% Inputs:
%   xy_obs   - n_obs x 2 matrix, where column 1 is time, column 2 is spatial coordinate
%   mCovF1   - covariance function handling
%   c1_t     - Parameters for temporal covariance function mCovF1
%   c1_s     - Parameters for spatial signal covariance function mCovF1
%   A_noise  - Amplitude for spatial noise covariance
%   t_master - Time of the master epoch for low variance assignment
%
% Output:
%   C_obs_obs_fin - n_obs x n_obs combined covariance matrix (signal + noise)
%   C_obs_obs     - n_obs x n_obs covariance matrix (signal)
%   C_vv_obs      - n_obs x n_obs covariance matrix (noise)

% Number of observations
n_obs = size(xy_obs, 1);

% 1. Calculate Distance Matrices (Vectorized)
% Create coordinate vectors
t_vec = xy_obs(:, 1);
s_vec = xy_obs(:, 2);

% Compute full grid of lags using implicit expansion or repmat
Tau_t = abs(t_vec - t_vec');
Tau_s = abs(s_vec - s_vec');

% 2. Evaluate Covariance Functions on Matrices
Ct_mat = mCovF1(c1_t, Tau_t);
Cs_mat = mCovF2(c1_s, Tau_s);

% 3. Combine (Separable Model)
% C(t,s) = C(t)*C(s) / sigma_t^2 to normalize amplitude
C_obs_obs = (Ct_mat .* Cs_mat) ./ c1_t(1);

% 4. Construct Noise Matrix (Diagonal)
% Create a vector of noise variances
diag_noise = repmat(A_noise, n_obs, 1);

% Reduce noise for master epoch
idx_master = (t_vec == t_master);
if any(idx_master)
    diag_noise(idx_master) = A_noise * 0.1;
end

% Create the sparse diagonal matrix (saves memory) or full matrix
Cvv_obs = diag(diag_noise);

% 5. Final Sum
C_obs_obs_fin = C_obs_obs + Cvv_obs;

% % 6. Eigenvalue check (Debug)
%     e = eig(C_obs_obs_fin);
%     e = sort(e, 'descend');
%     fprintf('Min eigenvalue   : %g\n', min(e));
%     fprintf('Max eigenvalue   : %g\n', max(e));
%     fprintf('Condition number : %.3e\n', max(e)/min(e));
%     fprintf('Negative eigenvalues: %d\n', sum(e < 0));
end