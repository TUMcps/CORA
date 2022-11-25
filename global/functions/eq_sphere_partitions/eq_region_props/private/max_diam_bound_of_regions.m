function diam_bound = max_diam_bound_of_regions(regions)
%MAX_DIAM_BOUND_OF_REGIONS The maximum diameter bound in an array of regions
%
% diam_bound = max_diam_bound_of_regions(regions);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function s2e changed name to sph2euc_dist
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(regions,1);
if dim == 1
    diam_bound = diam_bound_region(regions(:,:,1));
else
    colatitude = -inf*ones(dim-1,1);
    diam_bound = 0;
    N = size(regions,3);
    for region_n = 1:N
        top = regions(:,1,region_n);
        if norm(top(2:dim)-colatitude) ~= 0
            colatitude = top(2:dim);
            diam_bound = max(diam_bound,diam_bound_region(regions(:,:,region_n)));
        end
    end
end
%
% end function

function diam_bound = diam_bound_region(region)
% Calculate the per-region bound on the Euclidean diameter of a region.
%
% diam_bound = diam_bound_region(region)

tol = eps*2^5;
pseudo_region = pseudo_region_for_diam(region);
dim = size(pseudo_region,1);
top = pseudo_region(:,1);
bot = pseudo_region(:,2);
diam_bound = 0;
s = bot(dim)-top(dim);
e = sph2euc_dist(s);
if dim == 1
    diam_bound = e;
else
    max_sin = max(sin(top(dim)),sin(bot(dim)));
    if (top(dim) <= pi/2) && (bot(dim) >= pi/2)
        max_sin = 1;
    end
    if (abs(top(dim)) < tol) || (abs(pi-bot(dim)) < tol)
        diam_bound = 2*max_sin;
    else
        region_1 = [top(1:dim-1),bot(1:dim-1)];
        diam_bound_1 = max_sin*diam_bound_region(region_1);
        diam_bound = min(2,sqrt(e^2+diam_bound_1^2));
    end
end
%
% end function
