function coeff = calc_dist_coeff(dim,N,min_euclidean_dist)
%CALC_DIST_COEFF Coefficient of minimum distance
%
%Syntax
% coeff = calc_dist_coeff(dim,N,min_euclidean_dist);
%
%Description
% COEFF = CALC_DIST_COEFF(dim,N,MIN_EUCLIDEAN_DIST) sets COEFF to be the
% coefficient in the expression for the lower bound on the minimum distance of
% a minimum energy point set:
%
%    MIN_EUCLIDEAN_DIST >= COEFF N^(-1/dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The argument MIN_EUCLIDEAN_DIST must be an array of real nubers of the same array size as N.
% The result COEFF will be an array of the same size as N.
%
%Notes
% The expression for the lower bound on minimum distance of a minimum r^(-s)
% energy point set on S^dim was given by [RakSZ95] for s == 0 and dim = 2, 
% [Dahl78] for s == dim-1, [KuiSS04 Theorem 8] for dim-1 <= s < dim and
% [KuiS98 (1.12) p. 525] for s > dim.
%
%Examples
% > N=2:6
%  N =
%       2     3     4     5     6
%  
% > dist=eq_min_dist(2,N)
%  dist =
%      2.0000    1.4142    1.4142    1.4142    1.4142
%  
% > calc_dist_coeff(2,N,dist)
%  ans =
%      2.8284    2.4495    2.8284    3.1623    3.4641
%
%See also
% EQ_MIN_DIST, EQ_DIST_COEFF

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
error(nargchk(3,3,nargin));
%
% dim is the number of dimensions
% N is the number of regions
%
coeff =  min_euclidean_dist .* N.^(1/dim);
%
% end function
