function c_polar = polar_colat(dim, N)
%POLAR_COLAT The colatitude of the North polar (top) spherical cap
%
% Given dim and N, determine the colatitude of the North polar spherical cap.
%
% c_polar = polar_colat(dim, N);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

enough = N > 2;
c_polar(N == 1) = pi;
c_polar(N == 2) = pi/2;
c_polar(enough) = sradius_of_cap(dim,area_of_ideal_region(dim,N(enough)));
%
% end function
