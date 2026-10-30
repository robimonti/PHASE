function z=DSC_quadratic_example(x,y,s,t,p) %#ok<INUSD>
%DSC_QUADRATIC_EXAMPLE Rounded coefficients supplied in DSC.txt.
% Select lonlat coordinates, origin X=152.20890917, Y=-4.23823009.
% Model t is days since 2021-01-02; set the time offset if input starts later.
% This illustrative file does not assert the precision of the original fit.
c=polyval([2.5063e-9 -2.5302e-6 7.9440e-4 -9.1404e-2 36.816],t);
a=polyval([-1.7911e-7 2.9430e-4 -1.4059e-1 15.275 777.24],t);
b=polyval([2.5229e-7 -1.5958e-4 6.0671e-3 -4.291 2331.1],t);
d=polyval([6.2785e-6 -1.1970e-2 6.5309 272.05 -265470],t);
e=polyval([-1.1227e-5 1.4790e-2 -6.4103 884.75 63651],t);
f=polyval([4.8031e-6 -5.5241e-4 -1.9996 1222.1 -168600],t);
z=c+a.*x+b.*y+d.*x.*y+e.*x.^2+f.*y.^2;
end
