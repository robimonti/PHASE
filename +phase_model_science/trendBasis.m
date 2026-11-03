function B=trendBasis(times,reference,scale,degree)
%TRENDBASIS Polynomial displacements share the observations' reference epoch.
u=double(times(:))/scale;u0=reference/scale;
B=u.^(1:degree)-u0.^(1:degree);
end
