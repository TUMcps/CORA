function region = sphere_region(dim)
%SPHERE_REGION The sphere represented as a single region of an EQ partition
%
% An array of two points representing S^dim as a region.
%
% region = sphere_region(dim);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if dim == 1
    region = [0,2*pi];
else
    sphere_region_1 = sphere_region(dim-1);
    region = [[sphere_region_1(:,1); 0],[sphere_region_1(:,2); pi]] ;
end
%
% end function
