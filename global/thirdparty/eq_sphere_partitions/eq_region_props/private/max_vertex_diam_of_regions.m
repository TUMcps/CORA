function vertex_diam = max_vertex_diam_of_regions(regions)
%MAX_VERTEX_DIAM_OF_REGIONS The max vertex diameter in a cell array of regions
%
% vertex_diam = max_vertex_diam_of_regions(regions);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(regions,1);
if dim == 1
    vertex_diam = vertex_diam_region(regions(:,:,1));
else
    colatitude = -inf*ones(dim-1,1);
    vertex_diam = 0;
    N = size(regions,3);
    for region_n = 1:N
        top = regions(:,1,region_n);
        if norm(top(2:dim)-colatitude) ~= 0
            colatitude = top(2:dim);
            vertex_diam = max(vertex_diam,vertex_diam_region(regions(:,:,region_n)));
        end
    end
end
%
% end function

function diam = vertex_diam_region(region)
% Calculate the Euclidean diameter of the set of 2^dim vertices
% of the pseudo-region of a region.
%
% diam = vertex_diam_region(region);

expanded_region = expand_region_for_diam(region);
dim = size(expanded_region,1);
full = size(expanded_region,2);
half = floor(full/2);
top = expanded_region(:,1);
bot = expanded_region(:,full);
diam = 0;
if sin(top(dim)) > sin(bot(dim))
    for point_n_1 = 1:2:half
        for point_n_2 = point_n_1+1:2:full
            x1 = polar2cart(expanded_region(:,point_n_1));
            x2 = polar2cart(expanded_region(:,point_n_2));
            diam = max(diam,euclidean_dist(x1,x2));
        end
    end
else
    for point_n_1 = full:-2:half+1
        for point_n_2 = point_n_1-1:-2:1
            x1 = polar2cart(expanded_region(:,point_n_1));
            x2 = polar2cart(expanded_region(:,point_n_2));
            diam = max(diam,euclidean_dist(x1,x2));
        end
    end
end
%
% end function
