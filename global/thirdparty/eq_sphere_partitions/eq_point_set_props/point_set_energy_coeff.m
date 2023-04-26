function coeff = point_set_energy_coeff(points,s)
%POINT_SET_ENERGY_COEFF Coefficient in expansion of energy of a point set
%
%Syntax
% coeff = point_set_energy_coeff(points,s);
%
%Description
% COEFF = POINT_SET_ENERGY_COEFF(POINTS,s)  does the following:
% 1) finds the r^(-s) energy of the point set POINTS, and
% 2) sets COEFF to be the coefficient of the second term of the expansion of
%    the energy, having the same form as the expansion of E(dim,N,s), 
%    the minimum r^(-s) energy of a set of N points on S^dim.
%
% The argument POINTS must be an array of real numbers of size (dim+1 by N),
% where dim and N are positive integers. 
% Each column of POINTS represents a point in R^(dim+1).
%
% Specifically, for s not equal to 0, COEFF is the solution to
%
% energy == (SPHERE_INT_ENERGY(dim,s)/2) N^2 + COEFF N^(1+s/dim),
%
% and for s == 0 (the logarithmic potential), COEFF is the solution to
%
% energy == (SPHERE_INT_ENERGY(dim,0)/2) N^2 + COEFF N LOG(N).
%
% COEFF = POINT_SET_ENERGY_COEFF(POINTS) uses the default value dim-1 for s.
%
%Notes
% The value dim is the dimension of S^dim as a manifold. The point set POINTS
% is assumed to be a subset of R^(dim+1) but is not assumed to be a subset of
% S^dim.
%
% For details of the calculation of the energy coefficient, 
% see HELP CALC_ENERGY_COEFF.
%
%Examples
% > x
%  x =
%           0    0.0000   -0.0000    0.0000
%           0    1.0000   -1.0000         0
%      1.0000    0.0000    0.0000   -1.0000
%
% > point_set_energy_coeff(x)
%  ans =
%     -0.5214
%
%See also
% POINT_SET_ENERGY_DIST, CALC_ENERGY_COEFF, EQ_ENERGY_COEFF, EQ_ENERGY_DIST

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
error(nargchk(1,2,nargin));
%
% dim is the dimension of S^dim as a manifold.
%
dim = size(points,1)-1;
%
% N is the number of points in the point set.
%
N = size(points,2);
%
% The default value of s is dim-1.
%
if nargin < 2
    s = dim-1;
end

energy = point_set_energy_dist(points,s);
coeff = calc_energy_coeff(dim,N,s,energy);
%
% end function

