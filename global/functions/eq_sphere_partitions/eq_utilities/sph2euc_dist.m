function e = sph2euc_dist(s)
%SPHE2EUC_DIST Convert spherical to Euclidean distance
%
%Syntax
% e = sph2euc_dist(s);
%
%Description
% E = SPHE2EUC_DIST (S) converts the spherical distance S to Euclidean
% distance E, using a formula which is valid for the unit sphere in all
% dimensions.
%
% The argument S must be a real number or an array of real numbers.
% The result E will be an array of the same size as S.
%
%Note
% The argument S is assumed to satsify abs(S) <= pi.
%
%Examples
% > e=sph2euc_dist(pi)
% e =
%      2
%
% > e=sph2euc_dist(0:pi/4:pi)
% e =
%          0    0.7654    1.4142    1.8478    2.0000
%
%See also
% EUC2SPH_DIST, EUCLIDEAN_DIST, SPHERICAL_DIST

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Change name from s2e to sph2euc_dist
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

e = 2*(sin(s/2));
%
% end function
