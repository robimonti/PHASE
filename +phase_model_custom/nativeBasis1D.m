function B=nativeBasis1D(q,state)
% Original spline evaluation design from the automatic variance stage.
knots_x=state.knots_x; knots_y=state.knots_y;
n_params_spl=(numel(knots_x))*(numel(knots_y));
deltaX=knots_x(2)-knots_x(1); deltaY=knots_y(1)-knots_y(2);
phi_3=@(c) ((2-c).^3-4*(1-c).^3)/6; phi_4=@(c)(2-c).^3/6;
phi_33=@(a,b)phi_3(a).*phi_3(b); phi_34=@(a,b)phi_3(a).*phi_4(b);
phi_43=@(a,b)phi_4(a).*phi_3(b); phi_44=@(a,b)phi_4(a).*phi_4(b);
xy_est=q(:,[3,4]);
% 5.3.2) Build A_spline_est (prediction Grid)
x_est_vec = xy_est(:,1);
y_est_vec = xy_est(:,2);

n_points_est = length(x_est_vec);
max_nnz = n_points_est * 16;
I_idx = zeros(max_nnz, 1);
J_idx = zeros(max_nnz, 1);
V_val = zeros(max_nnz, 1);
count_nnz = 0;

for n = 1:n_points_est
    x = x_est_vec(n);
    y = y_est_vec(n);
    
    % Find knot interval for x
    i_x = find(x >= knots_x(1:end-1) & x < knots_x(2:end), 1);
    if isempty(i_x)
        if x < knots_x(1)
            i_x = 0;
        elseif x >= knots_x(end)
            i_x = length(knots_x) - 2;
        end
    else
        i_x = i_x - 1; 
    end
    
    % Find knot interval for y
    i_y = find(y <= knots_y(1:end-1) & y > knots_y(2:end), 1);
    if isempty(i_y)
        if y > knots_y(1)
            i_y = 0;
        elseif y <= knots_y(end)
            i_y = length(knots_y) - 2;
        end
    else
        i_y = i_y - 1; 
    end
    
    if (i_x >= -2 && i_x <= length(knots_x) - 1 && i_y >= -2 && i_y <= length(knots_y) - 1)
        % Compute local coordinates (FIXED)
        csi_x = (x - knots_x(i_x + 1)) / deltaX;
        csi_y = (knots_y(i_y + 1) - y) / deltaY;
        
        alpha = zeros(4, 4);
        alpha(1,1) = phi_44(1 + csi_x, 1 + csi_y); alpha(1,2) = phi_43(1 + csi_x, csi_y); alpha(1,3) = phi_43(1 + csi_x, 1 - csi_y); alpha(1,4) = phi_44(1 + csi_x, 2 - csi_y);
        alpha(2,1) = phi_34(csi_x, 1 + csi_y);     alpha(2,2) = phi_33(csi_x, csi_y);     alpha(2,3) = phi_33(csi_x, 1 - csi_y);     alpha(2,4) = phi_34(csi_x, 2 - csi_y);
        alpha(3,1) = phi_34(1 - csi_x, 1 + csi_y); alpha(3,2) = phi_33(1 - csi_x, csi_y); alpha(3,3) = phi_33(1 - csi_x, 1 - csi_y); alpha(3,4) = phi_34(1 - csi_x, 2 - csi_y);
        alpha(4,1) = phi_44(2 - csi_x, 1 + csi_y); alpha(4,2) = phi_43(2 - csi_x, csi_y); alpha(4,3) = phi_43(2 - csi_x, 1 - csi_y); alpha(4,4) = phi_44(2 - csi_x, 2 - csi_y);
        
        for k = -1:2
            for h = -1:2
                if (i_x + k >= 0 && i_x + k <= length(knots_x) - 1 && i_y + h >= 0 && i_y + h <= length(knots_y) - 1)
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

A_spline_est = sparse(I_idx(1:count_nnz), J_idx(1:count_nnz), V_val(1:count_nnz), n_points_est, n_params_spl);


B=A_spline_est;
end
