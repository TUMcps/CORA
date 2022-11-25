function n_regions = round_to_naturals(N,r_regions)
%ROUND_TO_NATURALS Round off a given list of numbers of regions
%
% Given N and r_regions, determine n_regions,
% a list of the natural number of regions in each collar and the polar caps.
% This list is as close as possible to r_regions, using rounding.
% The number of elements is n_collars+2.
% n_regions[1] is 1.
% n_regions[n_collars+2] is 1.
% The sum of n_regions is N.
%
% n_regions = round_to_naturals(N,r_regions);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_regions = r_regions;
discrepancy = 0;
for zone_n = 1:size(r_regions,2)
    n_regions(zone_n) = round(r_regions(zone_n)+discrepancy);
    discrepancy = discrepancy+r_regions(zone_n)-n_regions(zone_n);
end
%
% end function
