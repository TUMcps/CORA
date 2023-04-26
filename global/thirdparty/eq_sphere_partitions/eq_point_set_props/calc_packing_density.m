function density = calc_packing_density(dim,N,min_euclidean_dist)
%CALC_PACKING_DENSITY Density of packing given by minimum distance
%
%Syntax
% density = calc_packing_density(dim,N,min_euclidean_dist);
%
%Description
% DENSITY = CALC_PACKING_DENSITY(dim,N,MIN_EUCLIDEAN_DIST) sets DENSITY to
% be the density of a packing of S^dim by N equal spherical caps with centers
% having minimum Euclidean distance MIN_EUCLIDEAN_DIST.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The argument MIN_EUCLIDEAN_DIST must be an array of real nubers of the same array size as N.
% The result DENSITY will be an array of the same size as N.
%
%Notes
% The packing density is defined to be the sum of the areas of the spherical
% caps of the packing, divided by the area of the unit sphere S^dim.
%
% The spherical radius of the caps in the packing is half the minimum spherical
% distance between center points. The spherical radius for N == 1 is a special
% case. It is defined to be pi.
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
% > calc_packing_density(2,N,dist)
%  ans =
%      1.0000    0.4393    0.5858    0.7322    0.8787
%
%See also
% EQ_MIN_DIST, AREA_OF_CAP, AREA_OF_SPHERE, EQ_PACKING_DENSITY

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.10 $ $Date 2005-05-27 $
% Function e2s changed name to euc2sph_dist
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
s_cap = euc2sph_dist(min_euclidean_dist)/2;
%
% The spherical radius for N == 1 is a special case. It is pi.
%
s_cap(N == 1) = pi;
density = N.*area_of_cap(dim, s_cap)./area_of_sphere(dim);
%
% end function
