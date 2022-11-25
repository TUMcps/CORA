function coeff = eq_energy_coeff(dim,N,s,varargin)
%EQ_ENERGY_COEFF Coefficient in expansion of energy of an EQ point set
%
%Syntax
% coeff = eq_energy_coeff(dim,N,s,options);
%
%Description
% COEFF = EQ_ENERGY_COEFF(dim,N,s)  does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) finds the r^(-s) energy of the EQ point set, and
% 4) sets COEFF to be the coefficient of the second term of the expansion of
%    the energy, having the same form as the expansion of E(dim,N,s), 
%    the minimum r^(-s) energy of a set of N points on S^dim.
%
% Specifically, for s not equal to 0, COEFF is the solution to
%
% energy == (SPHERE_INT_ENERGY(dim,s)/2) N^2 + COEFF N^(1+s/dim),
%
% and for s == 0 (the logarithmic potential), COEFF is the solution to
%
% energy == (SPHERE_INT_ENERGY(dim,0)/2) N^2 + COEFF N LOG(N).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result COEFF will be an array of the same size as N.
%
% COEFF = EQ_ENERGY_COEFF(dim,N) uses the default value dim-1 for s.
%
% COEFF = EQ_ENERGY_COEFF(dim,N,s,'offset','extra') uses experimental extra offsets
% for S^2 and S^3 to try to minimize energy. For dim > 3, extra offsets are
% not used.
%
%Notes
% 1) The energy expansion is not valid for N == 1, and in particular,
%
% EQ_ENERGY_COEFF(dim,N,0) := 0.
%
% 2) For details of calculation of the energy coefficient, 
% see HELP CALC_ENERGY_COEFF
%
%Examples
% > coeff=eq_energy_coeff(2,10)
%  coeff =
%     -0.5461
%  
% > coeff=eq_energy_coeff(3,1:6)
%  coeff =
%     -0.5000   -0.5512   -0.5208   -0.5457   -0.5472   -0.5679
% > coeff=eq_energy_coeff(2,1:6,0)
%  coeff =
%           0   -0.2213   -0.1569   -0.2213   -0.2493   -0.2569
%
%See also
% PARTITION_OPTIONS, EQ_ENERGY_DIST, CALC_ENERGY_COEFF

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
error(nargchk(2,5,nargin));
%
% dim is the number of dimensions
% N is the number of regions
%
% The default value of s is dim-1.
%
if nargin < 3
    s = dim-1;
end
energy = eq_energy_dist(dim,N,s,varargin{:});
coeff = calc_energy_coeff(dim,N,s,energy);
%
% end function

