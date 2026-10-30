function B=nativeBasis2D(q,state)
% Original spline evaluation design from the automatic variance stage.
names=fieldnames(state);
for variableIndex=1:numel(names)
    eval([names{variableIndex} '=state.(names{variableIndex});']);
end
x_coords=q(:,1);y_coords=q(:,2);t_coords=q(:,4);
n_params=(num_x+2)*(num_y+2)*(num_t+2); n_ctrl_x=num_x+2;n_ctrl_y=num_y+2;n_ctrl_t=num_t+2;
deltaX=knots_x(2)-knots_x(1);deltaY=knots_y(2)-knots_y(1);deltaT=knots_t(2)-knots_t(1);
% 5.2.1) Rebuild the design matrix (A_obs) for the best model
% Needed to form the normal matrix (N) and get the covariance (C_cc)
n_obs_var = length(x_coords); 

% Triplets for A_obs_opt
max_nnz = n_obs_var * 64;
I_idx = zeros(max_nnz, 1);
J_idx = zeros(max_nnz, 1);
V_val = zeros(max_nnz, 1);
count_nnz = 0;

% Pre-compute inverse deltas
inv_dX = 1/deltaX;
inv_dY = 1/deltaY;
inv_dT = 1/deltaT;

for n = 1:n_obs_var
    x = x_coords(n);
    y = y_coords(n);
    t = t_coords(n);
    
    % Basis indices
    i_x = floor((x - knots_x(1)) * inv_dX);
    if i_x >= num_x, i_x = num_x - 1; end
    if i_x < 0, i_x = 0; end
    
    i_y = floor((y - knots_y(1)) * inv_dY);
    if i_y >= num_y, i_y = num_y - 1; end
    if i_y < 0, i_y = 0; end
    
    i_t = floor((t - knots_t(1)) * inv_dT);
    if i_t >= num_t, i_t = num_t - 1; end
    if i_t < 0, i_t = 0; end
    
    % Local coords
    csi_x = (x - knots_x(i_x + 1)) * inv_dX;
    csi_y = (y - knots_y(i_y + 1)) * inv_dY;
    csi_t = (t - knots_t(i_t + 1)) * inv_dT;
    
    % Basis functions
    bx = bspline_basis(i_x, csi_x);
    by = bspline_basis(i_y, csi_y);
    bt = bspline_basis(i_t, csi_t);
    
    % Fill Row a Triplette
    for l = -1:2
        idx_t = i_t + l + 1;
        for h = -1:2
            idx_y = i_y + h + 1;
            for k = -1:2
                idx_x = i_x + k + 1;
                if (idx_x >= 1 && idx_x <= n_ctrl_x && ...
                    idx_y >= 1 && idx_y <= n_ctrl_y && ...
                    idx_t >= 1 && idx_t <= n_ctrl_t)
                    
                    col = idx_y + (idx_x - 1) * n_ctrl_y + (idx_t - 1) * n_ctrl_y * n_ctrl_x;
                    
                    count_nnz = count_nnz + 1;
                    I_idx(count_nnz) = n;
                    J_idx(count_nnz) = col;
                    V_val(count_nnz) = bx(k+2) * by(h+2) * bt(l+2);
                end
            end
        end
    end
end

% Creation of the final A_obs_opt matrix
A_obs_opt = sparse(I_idx(1:count_nnz), J_idx(1:count_nnz), V_val(1:count_nnz), n_obs_var, n_params);


B=A_obs_opt;
end
