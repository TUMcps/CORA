function points_s = eq_point_set_polar(dim,N,varargin)
%EQ_POINT_SET_POLAR Center points of regions of an EQ partition
%
%Syntax
% points_s = eq_point_set_polar(dim,N,options);
%
%Description
% POINTS_S = EQ_POINT_SET_POLAR(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
% partition S^dim (the unit sphere in dim+1 dimensional space) into N regions
% of equal area and small diameter, and
% 2) sets POINTS_S to be an array of size (dim by N), containing the center
% points of each region. Each column of POINTS_S represents a point of S^dim,
% in spherical polar coordinates.
%
% The arguments dim and N must be positive integers.
%
% POINTS_S = EQ_POINT_SET_POLAR(dim,N,'offset','extra') uses experimental extra
% offsets for S^2 and S^3 to try to minimize energy. If dim > 3, extra offsets
% are not used.
%
% POINTS_S = EQ_POINT_SET_POLAR(dim,N,extra_offset) uses experimental extra
% offsets if extra_offset is true or non-zero.
%
%Notes
% Each region is defined as a product of intervals in spherical polar
% coordinates. The center point of a region is defined via the center points
% of each interval, with the exception of spherical caps and their descendants,
% where the center point is defined using the center of the spherical cap.
%
% For more details on options, see HELP PARTITION_OPTIONS.
%
%Examples
% > points_s = eq_point_set_polar(2,4)
% points_s =
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
% > size(points_s)
% ans =
%      2     4
%
%See also
% PARTITION_OPTIONS, EQ_POINT_SET, EQ_REGIONS

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Function changed name from x2s2 to cart2polar2
% Optimize running time:
%   use slice assignments
%   trade space for time by using a cache
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
error(nargchk(2,4,nargin));
%
% dim is the number of dimensions
% N is the number of regions
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
    fprintf('Usage: eq_point_set_polar(dim, N)\n');
    error('The arguments dim and N must be positive integers.');
end

if N == 1
    %
    % We have only one region, which must be the whole sphere.
    %
    points_s = zeros(dim,1);
    return;
end
%
% Start the partition of the sphere into N regions by partitioning
% to caps defined in the current dimension.
%
[a_cap, n_regions] = eq_caps(dim,N);
%
% a_cap is an increasing list of angles of the caps.
%
if dim == 1
    %
    % We have a circle and a_cap is an increasing list of angles of sectors,
    % with a_cap(k) being the cumulative arc length 2*pi/k.
    % The points are placed half way along each sector.
    %
    points_s = a_cap - pi/N;
    %
else
    %
    % We have a number of zones: two polar caps and a number of collars.
    % n_regions is the list of the number of regions in each zone.
    %
    n_collars = size(n_regions,2)-2;
    use_cache = dim >= 2;
    if use_cache
        cache_size = floor(n_collars/2);
        cache = cell(1,cache_size);
    end
    %
    % Start with the 'centre' point of the North polar cap.
    % This is the North pole.
    %
    points_s = zeros(dim,N);
    point_n = 2;
    %
    % Determine the 'centre' points for each collar.
    %
    if extra_offset && (dim == 3)
        R = eye(3);
    end
    if dim == 2
        offset = 0;
    end
    for collar_n = 1:n_collars
        %
        % a_top is the colatitude of the top of the current collar.
        %
        a_top = a_cap(collar_n);
        %
        % a_bot is the colatitude of the bottom of the current collar.
        %
        a_bot = a_cap(1+collar_n);
        %
        % n_in_collar is the number of regions in the current collar.
        %
        n_in_collar = n_regions(1+collar_n);
        %
        % The top and bottom of the collar are small (dim-1)-spheres,
        % which must be partitioned into n_in_collar regions.
        % Use eq_point_set_polar recursively to partition
        % the unit (dim-1)-sphere.
        % points_1 is the resulting list of points.
        %
        if use_cache
            twin_collar_n = n_collars-collar_n+1;
            if twin_collar_n <= cache_size && ...
                size(cache{twin_collar_n},2) == n_in_collar
                points_1 = cache{twin_collar_n};
            else
                points_1 = eq_point_set_polar(dim-1,n_in_collar,extra_offset);
                cache{collar_n} = points_1;
            end
        else
            points_1 = eq_point_set_polar(dim-1,n_in_collar,extra_offset);
        end
        %
        if extra_offset && (dim == 3) && (collar_n > 1)
            %
            % (Experimental)
            % Rotate 2-spheres to prevent alignment of north poles.
            %
            R = s2_offset(points_1)*R;
            points_1 = cart2polar2(R*polar2cart(points_1));
        end
        %
        % Given points_1, determine the 'centre' points for the collar.
        % Each point of points_1 is a 'centre' point on the (dim-1)-sphere.
        %
        % Angular 'centre' point;
        % The first angles are those of the current 'centre' point
        % of points_1, and the last angle in polar coordinates is the average of
        % the top and bottom angles of the collar,
        %
        a_point = (a_top+a_bot)/2;
        %
        point_1_n = 1:size(points_1,2);
        %
        if dim == 2
            %
            % The (dim-1)-sphere is a circle
            %
            points_s(1:dim-1,point_n+point_1_n-1) = mod(points_1(:,point_1_n)+2*pi*offset,2*pi);
            %
            % Given the number of sectors in the current collar and
            % in the next collar, calculate the next offset.
            % Accumulate the offset, and force it to be a number between 0 and 1.
            %
            offset = offset + circle_offset(n_in_collar,n_regions(2+collar_n),extra_offset);
            offset = offset - floor(offset);
        else
            points_s(1:dim-1,point_n+point_1_n-1) = points_1(:,point_1_n);
        end
        %
        points_s(dim, point_n+point_1_n-1) = a_point;
        point_n = point_n + size(points_1,2);
    end
    %
    % End with the 'centre' point of the bottom polar cap.
    %
    points_s(:,point_n) = zeros(dim,1);
    points_s(dim,point_n) = pi;
end
%
% end function
