function s = euc2sph_dist(e)
%EUC2SPH_DIST Convert Euclidean to spherical distance
%
%Syntax
% s = e2s(e);
%
%Description
% S = EUC2SPH_DIST(E) converts the Euclidean distance E to the spherical 
% distance S, using a formula which is valid for the unit sphere in all
% dimensions.
%
% The argument E must be a real number or an array of real numbers.
% The result S will be an array of the same size as E.
%
%Note
% The argument E is assumed to satsify abs(E) <= 2.
%
%Examples
% > s=euc2sph_dist(2)
%  s =
%      3.1416
%
% > s=euc2sph_dist(0:0.5:2)
%  s =
%           0    0.5054    1.0472    1.6961    3.1416
%
% > s=euc2sph_dist(-2)
%  s =
%     -3.1416
%
%See also
% SPH2EUC_DIST, EUCLIDEAN_DIST, SPHERICAL_DIST

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Change name from e2s to euc2sph_dist
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

s = 2*asin(e/2);
%
% end function
