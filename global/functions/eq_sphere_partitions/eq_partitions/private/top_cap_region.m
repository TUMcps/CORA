function region = top_cap_region(dim,a_cap)
%TOP_CAP_REGION North polar (top) cap region of EQ partition
%
% An array of two points representing the top cap of radius a_cap as a region.
%
% region = top_cap_region(dim,a_cap);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if dim == 1
    region = [0,a_cap];
else
    sphere_region_1 = sphere_region(dim-1);
    region = [[sphere_region_1(:,1); 0], [sphere_region_1(:,2); a_cap]];
end
% end function
