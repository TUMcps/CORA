function r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%IDEAL_REGION_LIST The ideal real number of regions in each zone
%
% List the ideal real number of regions in each collar, plus the polar caps.
%
% Given dim, N, c_polar and n_collars, determine r_regions,
% a list of the ideal real number of regions in each collar,
% plus the polar caps.
% The number of elements is n_collars+2.
% r_regions[1] is 1.
% r_regions[n_collars+2] is 1.
% The sum of r_regions is N.
%
% r_regions = ideal_region_list(dim,N,c_polar,n_collars);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

r_regions = zeros(1,2+n_collars);
r_regions(1) = 1;
if n_collars > 0
    %
    % Based on n_collars and c_polar, determine a_fitting,
    % the collar angle such that n_collars collars fit between the polar caps.
    %
    a_fitting = (pi-2*c_polar)/n_collars;
    ideal_region_area = area_of_ideal_region(dim,N);    
    for collar_n = 1:n_collars
        ideal_collar_area = area_of_collar(dim, c_polar+(collar_n-1)*a_fitting, c_polar+collar_n*a_fitting);
        r_regions(1+collar_n) = ideal_collar_area / ideal_region_area;
    end
end
r_regions(2+n_collars) = 1;

% end function
