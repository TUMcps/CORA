function region = bot_cap_region(dim,a_cap)
%BOT_CAP_REGION South polar (bottom) cap region of EQ partition
%
% An array of two points representing the bottom cap of radius a_cap as a region.
%
% region = bot_cap_region(dim,a_cap);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if dim == 1
    region = [2*pi-a_cap,2*pi];
else
    sphere_region_1 = sphere_region(dim-1);
    region = [[sphere_region_1(:,1); pi-a_cap],[sphere_region_1(:,2); pi]];
end
%
% end function
