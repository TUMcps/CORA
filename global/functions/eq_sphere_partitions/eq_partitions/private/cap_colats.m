function c_caps = cap_colats(dim,N,c_polar,n_regions)
%CAP_COLATS Colatitudes of spherical caps enclosing cumulative sum of regions
%
% Given dim, N, c_polar and n_regions, determine c_caps,
% an increasing list of colatitudes of spherical caps which enclose the same area
% as that given by the cumulative sum of regions.
% The number of elements is n_collars+2.
% c_caps[1] is c_polar.
% c_caps[n_collars+1] is Pi-c_polar.
% c_caps[n_collars+2] is Pi.
%
% c_caps = cap_colats(dim,N,c_polar,n_regions);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

c_caps = zeros(size(n_regions));
c_caps(1) = c_polar;
ideal_region_area = area_of_ideal_region(dim,N);
n_collars = size(n_regions,2)-2;
subtotal_n_regions = 1;
for collar_n = 1:n_collars
    subtotal_n_regions = subtotal_n_regions+n_regions(1+collar_n);
    c_caps(collar_n+1) =sradius_of_cap(dim,subtotal_n_regions*ideal_region_area);
end
c_caps(1+n_collars+1) = pi;

% end function
