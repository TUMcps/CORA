function density = point_set_packing_density(points)
%POINT_SET_PACKING_DENSITY  Density of packing given by minimum distance of a point set
%
%Syntax
% density = point_set_packing_density(points);
%
%Description
% DENSITY = POINT_SET_PACKING_DENSITY(POINTS) does the following:
% 1) finds the minimum Euclidean distance between points of the point set
%    POINTS, and
% 2) sets DENSITY to be the density of a packing of S^dim by N equal
% spherical caps with this minimum distance.
%
% The argument POINTS must be an array of real numbers of size (dim+1 by N),
% where dim and N are positive integers. 
% Each column of POINTS must represents a point of S^dim in Cartesian 
% coordinates.
%
%Notes
% Because packing density is defined using spherical caps, it well defined only
% for points on S^dim. Therefore POINTS must represent a subset of S^dim.
%
% For more details on the calculation of packing density, 
% see HELP CALC_PACKING_DENSITY.
%
%Examples 
% > x
%  x =
%           0    0.0000   -0.0000    0.0000
%           0    1.0000   -1.0000         0
%      1.0000    0.0000    0.0000   -1.0000
%  
% > point_set_packing_density(x)
%  ans =
%      0.5858
%
%See also
% CALC_PACKING_DENSITY, EQ_MIN_DIST, AREA_OF_CAP, AREA_OF_SPHERE

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
error(nargchk(1,1,nargin));
%
% Check that points lie on the unit sphere.
%
tol = eps * 2^5;
radius = sqrt(sum(points.*points));
if (min(radius) < 1-tol) || (max(radius) > 1+tol)
    error('Point set must be a subset of the unit sphere');
end    
%
% dim is the dimension of S^dim as a manifold.
%
dim = size(points,1)-1;
%
% N is the number of points in the point set.
%
N = size(points,2);
%
min_euclidean_dist = min(2,point_set_min_dist(points));
density = calc_packing_density(dim,N,min_euclidean_dist);
%
% end function
