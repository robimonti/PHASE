function val = bspline_basis(~, csi)

% bspline_basis Returns 4 values for the basis functions active at knot i
% corresponding to N_{i-1}, N_i, N_{i+1}, N_{i+2}
    
    % Using standard cardinal form for u in [0,1]:
    t = csi;
    val = zeros(4,1);
    
    % These are standard Uniform Cubic B-Spline basis functions
    % B_0 (for node i-1): (1-t)^3/6
    % B_1 (for node i)  : (3t^3 - 6t^2 + 4)/6
    % B_2 (for node i+1): (-3t^3 + 3t^2 + 3t + 1)/6
    % B_3 (for node i+2): t^3/6
    val(1) = (1 - t)^3 / 6;
    val(2) = (3*t^3 - 6*t^2 + 4) / 6;
    val(3) = (-3*t^3 + 3*t^2 + 3*t + 1) / 6;
    val(4) = t^3 / 6;
end