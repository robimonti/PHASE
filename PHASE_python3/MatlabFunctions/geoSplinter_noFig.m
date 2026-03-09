function [data, raster, par, sigmaPar, nodes_coord] = geoSplinter_noFig (fileName, method)

% SINTAX:
%   [data, raster, par, sigmaPar] = geoSplinter (fileName, method);
%
% INPUT:
%   fileName = input file name (without extension)
%   method   = interpolation method:
%              - 'lin' = linear splines    (1D)
%              - 'cub' = cubic splines     (1D)
%              - 'bil' = bilinear splines  (2D)
%              - 'bic' = bicubic splines   (2D)
%
% OUTPUT:
%   data     = observation points, observed and estimated values
%   raster   = estimated values on a regular grid
%   par      = spline parameters
%   sigmaPar = standard deviations of the spline parameters
%
% DESCRIPTION:
%   Representation of the interpolated data.
%
% SEE ALSO   ...
%
% Mirko Reguzzoni
% OGS - Politecnico di Milano
% 22-01-2005

%-----------------------------------------------------------------------

if ((method == 'lin') | (method == 'cub'))
	dim = 1;
else % if ((method == 'bil') | (method == 'bic'))
	dim = 2;
end

%-----------------------------------------------------------------------

data     = load ([fileName, '.out.txt']);
raster   = load ([fileName, '.grd.txt']);
par      = load ([fileName, '.par.txt']);
sigmaPar = load ([fileName, '.std.txt']);

%-----------------------------------------------------------------------

fid = fopen ([fileName, '.hdr.txt'], 'r');

str  = fscanf (fid, '%s %s', 2);
xMin = fscanf (fid, '%f', 1);

str  = fscanf (fid, '%s %s', 2);
xMax = fscanf (fid, '%f', 1);

if (dim == 2)

	str  = fscanf (fid, '%s %s', 2);
	yMin = fscanf (fid, '%f', 1);

	str  = fscanf (fid, '%s %s', 2);
	yMax = fscanf (fid, '%f', 1);

end

str    = fscanf (fid, '%s %s', 2);
deltaX = fscanf (fid, '%f', 1);

xNum = round((xMax - xMin) / deltaX) + 1;
deltaX = (xMax - xMin) / (xNum - 1);

if (dim == 2)

	str    = fscanf (fid, '%s %s', 2);
	deltaY = fscanf (fid, '%f', 1);

	yNum = round((yMax - yMin) / deltaY) + 1;
	deltaY = (yMax - yMin) / (yNum - 1);

end

% str    = fscanf (fid, '%s %s', 2);
% digits = fscanf (fid, '%f', 1);

fclose (fid);



%-----------------------------------------------------------------------

if (dim == 1)

   if (method == 'lin')

       nodes_coord = xMin:deltaX:xMax;

   else % if (method == 'cub')
		xx = xMin : deltaX : xMax;
		nodes_coord = xx;
   end

else % if (dim == 2)

   xx = xMin:deltaX:xMax;
   yy = yMin:deltaY:yMax;

   x_grid = repmat(xx, length(yy), 1);
   y_grid = repmat(fliplr(yy)', 1, length(xx));

   nodes_coord(:,:,1) = x_grid;
   nodes_coord(:,:,2) = y_grid;
   
end


% -------------------------------------------------------------------------
end