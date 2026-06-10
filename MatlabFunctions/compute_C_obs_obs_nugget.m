function [C_obs_obs_fin, C_obs_obs, Cvv_obs] = compute_C_obs_obs_nugget(xy_obs, mCovF1, c1_t, mCovF2, c1_s, A_noise, t_master)

% compute_C_obs_obs_nugget Computes the combined covariance matrix C_obs_obs for 1D Space + Time data.
%
% Inputs:
%   xy_obs   - n_obs x 2 matrix: [Time, Spatial coordinate]
%   mCovF1   - Function handle for temporal covariance
%   c1_t     - Parameters for temporal covariance
%   mCovF2   - Function handle for spatial covariance
%   c1_s     - Parameters for spatial covariance
%   A_noise  - Amplitude for spatial noise (Nugget)
%   t_master - Master epoch time (for variance reduction)

% Number of observations
n_obs = size(xy_obs, 1);

% 1. Extract Vectors
t_vec = xy_obs(:, 1);
s_vec = xy_obs(:, 2);

% --- OPTIMIZATION 1: FORCE SINGLE PRECISION FOR RAM ---
C_obs_obs = zeros(n_obs, n_obs, 'single');

% --- OPTIMIZATION 2: CHUNKING ---
blockSize = 2000; 

for i = 1:blockSize:n_obs
    idx_end = min(i + blockSize - 1, n_obs);
    row_idx = i:idx_end;
    
    % Vectors limited to current block
    t_block = t_vec(row_idx);
    s_block = s_vec(row_idx);
    
    % Lags calculations for the block (blockSize x n_obs) -> Takes very little RAM
    Tau_t_block = abs(t_block - t_vec');
    Tau_s_block = abs(s_block - s_vec');
    
    % Evaluate Covariance Functions on block matrices
    Ct_mat_block = mCovF1(c1_t, Tau_t_block);
    Cs_mat_block = mCovF2(c1_s, Tau_s_block);
    
    % Insert data in the pre-allocated single matrix
    if c1_t(1) ~= 0
        C_obs_obs(row_idx, :) = single((Ct_mat_block .* Cs_mat_block) ./ c1_t(1));
    else
        C_obs_obs(row_idx, :) = 0;
    end
end

if c1_t(1) == 0
    warning('Temporal variance is zero. Returning zero covariance.');
end

% 4. Construct Noise Matrix (Diagonal)
diag_noise = repmat(A_noise, n_obs, 1);

% Reduce noise for master epoch
idx_master = (t_vec == t_master);
if any(idx_master)
    diag_noise(idx_master) = A_noise * 0.1;
end

% --- OPTIMIZATION 3: CRITICAL FIX — SPARSE MATRIX ---
% CRUCIAL: Changed from dense diag() to sparse spdiags(). 
Cvv_obs = spdiags(single(diag_noise), 0, n_obs, n_obs);

% 5. Final Sum (Performed in-place to save memory)
C_obs_obs_fin = C_obs_obs;
idx_diag = 1:(n_obs+1):(n_obs^2);
C_obs_obs_fin(idx_diag) = C_obs_obs_fin(idx_diag) + single(diag_noise)';

end