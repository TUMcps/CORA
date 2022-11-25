function coeff = point_set_dist_coeff(points)
%POINT_SET_DIST_COEFF Coefficient of minimum distance of a point set
%
%Syntax
% coeff = point_set_dist_coeff(points);
%
%Description
% COEFF = POINT_SET_DIST_COEFF(POINTS) does the following:
% 1) finds the minimum Euclidean distance between points of the point set
%    POINTS, which should be a subset of the unit sphere S^dim, and
% 2) sets COEFF to be the coefficient in the expression for the lower bound on
%    the minimum distance of a minimum energy point set:
%
%    MIN_EUCLIDEAN_DIST >= COEFF N^(-1/dim),
%
%    where N is the number of points in POINTS.
%
% The argument POINTS must be an array of real numbers of size (dim+1 by N),
% where dim and N are positive integers. 
% Each column of POINTS represents a point in R^(dim+1).
% It is assumed that point set POINTS is a subset of the unit sphere S^dim,
% but this is not checked.
%
%Notes
% Fore more details on the calculation of the coefficient of the minimum
% distance, see HELP CALC_DIST-COEFF.
%
%Examples
% > x
%  x =
%           0    0.0000   -0.0000    0.0000
%           0    1.0000   -1.0000         0
%      1.0000    0.0000    0.0000   -1.0000
%  
% > point_set_dist_coeff(x)
%  ans =
%      2.8284
%
%See also
% POINT_SET_MIN_DIST, CALC_DIST_COEFF, EQ_DIST_COEFF, EQ_MIN_DIST

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
% dim is the dimension of S^dim as a manifold.
%
dim = size(points,1)-1;
%
% N is the number of points in the point set.
%
N = size(points,2);
%
min_euclidean_dist = point_set_min_dist(points);
coeff = calc_dist_coeff(dim,N,min_euclidean_dist);
%
% end function
