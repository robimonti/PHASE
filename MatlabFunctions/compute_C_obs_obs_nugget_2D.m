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

% Pre-allocation
C_obs_obs = zeros(n_obs, n_obs, 'single');

% 2. Calculate Distance Matrices (Chunking)
blockSize = 2000; 

for i = 1:blockSize:n_obs
    idx_end = min(i + blockSize - 1, n_obs);
    row_idx = i:idx_end;
    
    % Vectors limited to current block
    t_block = t_vec(row_idx);
    x_block = x_vec(row_idx);
    y_block = y_vec(row_idx);
    
    % Asymmetric block matrices (blockSize x n_obs)
    Tau_t_block = abs(t_block - t_vec');
    Diff_X_block = x_block - x_vec';
    Diff_Y_block = y_block - y_vec';
    Tau_s_block  = sqrt(Diff_X_block.^2 + Diff_Y_block.^2);
    
    % Block covariance
    Ct_mat_block = mCovF1(c1_t, Tau_t_block);
    Cs_mat_block = mCovF2(c1_s, Tau_s_block);
    
    % Insert data in pre-allocated matrix
    if c1_t(1) ~= 0
        C_obs_obs(row_idx, :) = single((Ct_mat_block .* Cs_mat_block) ./ c1_t(1));
    else
        C_obs_obs(row_idx, :) = 0;
    end
end

if c1_t(1) == 0
    warning('Temporal variance is zero. Returning zero covariance.');
end
    
% 5. Construct Noise Matrix
diag_noise = repmat(A_noise, n_obs, 1);

% Reduce noise for master epoch (constraint enforcement)
% Using a small tolerance for float comparison
idx_master = abs(t_vec - t_master) < 1e-5;
if any(idx_master)
    diag_noise(idx_master) = A_noise * 0.1;
end

% Create diagonal matrix
Cvv_obs = spdiags(diag_noise, 0, n_obs, n_obs);

% 6. Final Sum
C_obs_obs_fin = C_obs_obs;
idx_diag = 1:(n_obs+1):(n_obs^2);
C_obs_obs_fin(idx_diag) = C_obs_obs_fin(idx_diag) + diag_noise';

% 7. Eigenvalue check (Optional Debugging)
% if n_obs <= 2000
%     e = eig(C_obs_obs_fin);
%     if min(e) < 0
%         fprintf('Warning: Covariance matrix is not positive definite.\n');
%         fprintf('Min eigenvalue   : %g\n', min(e));
%     end
% end
end