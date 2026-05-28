function B = compute_B_spline_row(t, knots, deltaX, n_params)
    % Compute B-spline design matrix row for a given time t
    % t: time point (e.g., master_t_rel)
    % knots: knot vector from geoSplinter
    % deltaX: knot spacing
    % n_params: number of parameters (n_spl + 2 for cubic splines)

    phi_3 = @(csi) (csi >= 0 & csi < 1) .* ((2 - csi).^3 - 4 * (1 - csi).^3) / 6 + ...
               (csi >= 1 & csi < 2) .* ((2 - csi).^3) / 6;
    phi_4 = @(csi) (csi >= 0 & csi < 2) .* ((2 - csi).^3) / 6;


    B = zeros(1, n_params);
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
    alpha(1) = phi_4(csi_i + 1); % B_{i_i-1}
    alpha(2) = phi_3(csi_i);     % B_{i_i}
    alpha(3) = phi_3(1 - csi_i); % B_{i_i+1}
    alpha(4) = phi_4(2 - csi_i); % B_{i_i+2}
    for s = -1:2
        j = i_i + s + 1;
        if j >= 1 && j <= n_params
            B(j) = alpha(s + 2);
        end
    end
end