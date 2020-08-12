function coeff = eq_dist_coeff(dim,N,varargin)
%EQ_DIST_COEFF Coefficient of minimum distance of an EQ point set
%
%Syntax
% coeff = eq_dist_coeff(dim,N,options);
%
%Description
% COEFF = EQ_DIST_COEFF(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) finds the minimum Euclidean distance between points of the EQ point set,
% 4) sets COEFF to be the coefficient in the expression for the lower bound on
%    the minimum distance of a minimum energy point set:
%
%    DIST >= COEFF N^(-1/dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result COEFF will be an array of the same size as N.
%
% COEFF = EQ_DIST_COEFF(dim,N,'offset','extra'), for dim == 2 or dim == 3, uses
% experimental extra rotation offsets to try to maximize the minimum distance. 
% For dim > 3, extra offsets are not used.
%
%Notes
% The expression for the lower bound on minimum distance of a minimum r^(-s)
% energy point set on S^dim was given by [RakSZ95] for s == 0 and dim = 2, 
% [Dahl78] for s == dim-1, [KuiSS04 Theorem 8] for dim-1 <= s < dim and
% [KuiS98 (1.12) p. 525] for s > dim.
%
% Ideally eq_dist_coeff(dim,N) should tend to area_of_sphere(dim)^(1/dim) as 
% N goes to infinity.
%
%Examples
% > coeff=eq_dist_coeff(2,10)
%  coeff =
%      3.3250
%  
% > coeff=eq_dist_coeff(3,1:6)
%  coeff =
%      2.0000    2.5198    2.0396    2.2449    2.4183    2.5698
%
%See also
% PARTITION_OPTIONS, EQ_MIN_DIST

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
error(nargchk(2,4,nargin));
%
% dim is the number of dimensions
% N is the number of regions
%
dist = eq_min_dist(dim,N,varargin{:});
coeff =  dist .* N.^(1/dim);
%
% end function
