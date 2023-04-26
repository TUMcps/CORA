function [regions,dim_1_rot] = eq_regions(dim,N,varargin)
%EQ_REGIONS Recursive zonal equal area (EQ) partition of sphere
%
%Syntax
% [regions,dim_1_rot] = eq_regions(dim,N,options);
%
%Description
% REGIONS = EQ_REGIONS(dim,N) uses the recursive zonal equal area sphere
% partitioning algorithm to partition S^dim (the unit sphere in dim+1
% dimensional space) into N regions of equal area and small diameter.
%
% The arguments dim and N must be positive integers.
%
% The result REGIONS is a (dim by 2 by N) array, representing the regions
% of S^dim. Each element represents a pair of vertex points in spherical polar
% coordinates.
%
% Each region is defined as a product of intervals in spherical polar
% coordinates. The pair of vertex points regions(:,1,n) and regions(:,2,n) give
% the lower and upper limits of each interval.
%
% REGIONS = EQ_REGIONS(dim,N,'offset','extra') uses experimental extra 
% offsets for S^2 and S^3 to try to minimize energy. If dim > 3, extra offsets 
% are not used. 
%
% REGIONS = EQ_REGIONS(dim,N,extra_offset) uses experimental extra offsets 
% if extra_offset is true or non-zero.
%
% [REGIONS,DIM_1_ROT] = EQ_REGIONS(dim,N) also returns DIM_1_ROT, a cell 
% array containing N rotation matrices, one per region, each of size dim by dim.
% These describe the R^dim rotation needed to place the region in its final 
% position.
%
% [REGIONS,DIM_1_ROT] = EQ_REGIONS(dim,N,'offset','extra') partitions S^dim 
% into N regions, using extra offsets, and also returning DIM_1_ROT, as above.
%
%Notes
% The output argument DIM_1_ROT is the only way to track the effect of the extra
% offset when dim == 3, because the R^3 rotation means that the boundary of a 
% region generally no longer coincides with hyperplanes of colatitude and 
% longitude. The function ILLUSTRATE_S3_PARTITION uses DIM_1_ROT.
%
% For more details on options, see HELP PARTITION_OPTIONS.
%
%Examples
% > regions = eq_regions(2,4)
% regions(:,:,1) =
%          0    6.2832
%          0    1.0472
% regions(:,:,2) =
%          0    3.1416
%     1.0472    2.0944
% regions(:,:,3) =
%     3.1416    6.2832
%     1.0472    2.0944
% regions(:,:,4) =
%          0    6.2832
%     2.0944    3.1416
%
% > size(regions)
% ans =
%      2     2     4
%
%See also
% PARTITION_OPTIONS, EQ_POINT_SET, EQ_POINT_SET_POLAR, PROJECT_S3_PARTITION

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Fix bug in assignment of dim_1_rot
% Documentation files renamed
% $Revision 1.02 $ $Date 2005-04-16 $
% Optimize running time:
%   move 'if nargout' blocks, refactor slice assignments
%   trade space for time by using a cache
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
error(nargchk(2,4,nargin));
error(nargoutchk(0,2,nargout));
%
% If the option 'offset' is 'extra', then use experimental extra offsets
% for S^2 and S^3 regions to try to minimize energy.
% The default is not to use extra offsets.
%
pdefault.extra_offset =  false;
if nargin < 3
    extra_offset = pdefault.extra_offset;
else
    popt = partition_options(pdefault, varargin{:});
    extra_offset = popt.extra_offset;
end
%
% Extra offset does not currently work for dim > 3,
% so quietly ignore this option in this case.
% Note that this also affects recursive calls to lower dimensions.
%
if dim > 3
    extra_offset = false;
end
%
% Check that dim and N are positive integers.
%
if ~( isnumeric(dim) && (dim >= 1) && (dim == floor(dim)) ) || ...
   ~( isnumeric(N) && (N >= 1) && (N == floor(N)) )
    fprintf('Usage: eq_regions(dim, N)\n');
    error('The arguments dim and N must be positive integers.');
end

if nargout > 1
    dim_1_rot = cell(1,N);
end

if N == 1
    %
    % We have only one region, which must be the whole sphere.
    %
    regions = zeros(dim,2,1);
    regions(:,:,1) = sphere_region(dim);
    if nargout > 1
        dim_1_rot{1} = eye(dim);
    end
    return;
end
%
% Start the partition of the sphere into N regions by partitioning
% to caps defined in the current dimension.
%
[s_cap, n_regions] = eq_caps(dim,N);
%
% s_cap is an increasing list of colatitudes of the caps.
%
if dim == 1
    %
    % We have a circle and s_cap is an increasing list of angles of sectors.
    %
    if nargout > 1
        R = eye(dim);
        for region_n = 1:N
            dim_1_rot{region_n} = R;
        end
    end
    %
    % Return a list of pairs of sector angles.
    %
    regions = zeros(dim,2,N);
    regions(:,1,2:N) = s_cap(1:N-1);
    regions(:,2,:)   = s_cap;
    %
else
    %
    % We have a number of zones: two polar caps and a number of collars.
    % n_regions is the list of the number of regions in each zone.
    %
    n_collars = size(n_regions,2)-2;
    use_cache = dim > 2;
    if use_cache
        cache_size = floor(n_collars/2);
        cache = cell(1,cache_size);
    end
    %
    % Start with the top cap.
    %
    regions = zeros(dim,2,N);
    regions(:,:,1) = top_cap_region(dim,s_cap(1));
    region_n = 1;
    %
    % Determine the dim-regions for each collar.
    %
    if (nargout > 1) || (extra_offset && (dim == 3))
        R = eye(dim);
    end
    if nargout > 1
        dim_1_rot{1} = R;
    end
    if dim == 2
        offset = 0;
    end
    for collar_n = 1:n_collars
        %
        % c_top is the colatitude of the top of the current collar.
        %
        c_top = s_cap(collar_n);
        %
        % c_bot is the colatitude of the bottom of the current collar.
        %
        c_bot = s_cap(1+collar_n);
        %
        % n_in_collar is the number of regions in the current collar.
        %
        n_in_collar = n_regions(1+collar_n);
        %
        % The top and bottom of the collar are small (dim-1)-spheres,
        % which must be partitioned into n_in_collar regions.
        % Use eq_regions recursively to partition
        % the unit (dim-1)-sphere.
        % regions_1 is the resulting list of (dim-1)-region pairs.
        %
        if use_cache
            twin_collar_n = n_collars-collar_n+1;
            if twin_collar_n <= cache_size && ...
                size(cache{twin_collar_n},3) == n_in_collar
                regions_1 = cache{twin_collar_n};
            else
                regions_1 = eq_regions(dim-1,n_in_collar,extra_offset);
                cache{collar_n} = regions_1;
            end
        else
            regions_1 = eq_regions(dim-1,n_in_collar,extra_offset);
        end
        %
        if extra_offset && (dim == 3) && (collar_n > 1)
            %
            % (Experimental)
            % Rotate 2-spheres to prevent alignment of north poles.
            %
            R = s2_offset(centres_of_regions(regions_1))*R;
        end
        %
        % Given regions_1, determine the dim-regions for the collar.
        % Each element of regions_1 is a (dim-1)-region pair for
        % the (dim-1)-sphere.
        %
        if nargout > 1
            for region_1_n = 1:size(regions_1,3)
                dim_1_rot{region_n+region_1_n} = R;
            end
        end
        if dim == 2
            %
            % The (dim-1)-sphere is a circle
            % Offset each sector angle by an amount which accumulates over
            % each collar.
            %
            for region_1_n = 1:size(regions_1,3)
                %
                % Top of 2-region;
                % The first angle is the longitude of the top of
                % the current sector of regions_1, and
                % the second angle is the top colatitude of the collar.
                %
                r_top = [mod(regions_1(1,1,region_1_n)+2*pi*offset,2*pi); c_top];
                %
                % Bottom of 2-region;
                % The first angle is the longitude of the bottom of
                % the current sector of regions_1, and
                % the second angle is the bottom colatitude of the collar.
                %
                r_bot = [mod(regions_1(1,2,region_1_n)+2*pi*offset,2*pi); c_bot];
                if r_bot(1) < r_top(1)
                   r_bot(1) = r_bot(1) + 2*pi;
                end
                region_n = region_n+1;
                regions(:,:,region_n) = [r_top,r_bot];
            end
            %
            % Given the number of sectors in the current collar and
            % in the next collar, calculate the next offset.
            % Accumulate the offset, and force it to be a number between 0 and 1.
            %
            offset = offset + circle_offset(n_in_collar,n_regions(2+collar_n),extra_offset);
            offset = offset - floor(offset);
        else
            for region_1_n = 1:size(regions_1,3)
                %
                region_n = region_n+1;
                %
                % Dim-region;
                % The first angles are those of the current (dim-1) region of regions_1.
                %
                regions(1:dim-1,:,region_n) = regions_1(:,:,region_1_n);
                %
                % The last angles are the top and bottom colatitudes of the collar.
                %
                regions(dim,:,region_n) = [c_top,c_bot];
            end
        end
    end
    %
    % End with the bottom cap.
    %
    regions(:,:,N) = bot_cap_region(dim,s_cap(1));
    if nargout > 1
        dim_1_rot{N} = eye(dim);
    end
end
%
% end function
