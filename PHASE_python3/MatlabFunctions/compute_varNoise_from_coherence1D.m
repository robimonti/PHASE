function varNoise = compute_varNoise_from_coherence1D(coh_values, lambda_sar, num_looks, lambda_fallback)

% compute_varNoise_from_coherence1D
%
% Inputs:
%   coh_values: [N x M] matrix of coherence values per PS and epoch
%   lambda_sar: Radar wavelength (m, e.g., 0.05546576 for Sentinel-1)
%   num_looks: Number of looks (e.g., 4)
%   lambda_fallback: Fallback variance value in mm^2
%
% Output:
%   varNoise: [N x 1] vector of noise variance (mm^2) per PS

N_ps = size(coh_values, 1);
varNoise = zeros(N_ps, 1);

for i = 1:N_ps
    coh_valid = coh_values(i, :);
    coh_valid = coh_valid(~isnan(coh_valid));
    if ~isempty(coh_valid)
        coh_mean = median(coh_valid, 'omitnan');
        if coh_mean > 0.1
            sigma_phi2 = (1 - coh_mean^2) / (2 * coh_mean^2 * num_looks);
            varNoise(i) = ((lambda_sar / (4 * pi))^2) * sigma_phi2 * 1e6; % mm^2
        else
            varNoise(i) = lambda_fallback; % Fallback
        end
    else
        varNoise(i) = lambda_fallback; % Fallback
    end
end
end