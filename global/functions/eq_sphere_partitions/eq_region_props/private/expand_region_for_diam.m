function expanded_region = expand_region_for_diam(region)
%EXPAND_REGION_FOR_DIAM The set of 2^d vertices of a region
%
% Expand a region from the 2 vertex definition to the set of 2^dim vertices
% of the pseudo-region of a region, so that the Euclidean diameter of a region 
% is approximated by the diameter of this set.
%
% expanded_region = expand_region_for_diam(region);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(region,1);
if dim > 1
    s_top = region(dim,1);
    s_bot = region(dim,2);
    region_1 = expand_region_for_diam(region(1:dim-1,:));
    expanded_region = [append(region_1, s_top), append(region_1, s_bot)];
else
    expanded_region = pseudo_region_for_diam(region);
end
%
% end function

function result = append(matrix,value)
% Append a coordinate value to each column of a matrix.
%
% result = append(matrix,value);

result = [matrix; ones(1,size(matrix,2))*value];
%
% end function
