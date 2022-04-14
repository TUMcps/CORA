function [total_error, max_error] = eq_area_error(dim,N)
%EQ_AREA_ERROR Total area error and max area error per region of an EQ partition
%
%Syntax
% [total_error, max_error] = eq_area_error(dim,N)
%
%Description
% [TOTAL_ERROR, MAX_ERROR] = EQ_AREA_ERROR(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) sets TOTAL_ERROR to be the absolute difference between the total area of
%    all regions of the partition, and the area of S^dim, and
% 3) sets MAX_ERROR to be the maximum absolute difference between the area of 
%    any region of the partition, and the ideal area of a region as given by
%    AREA_OF_IDEAL_REGION(dim,N), which is 1/N times the area of S^dim.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The results TOTAL_ERROR and MAX_ERROR will be arrays of the same size as N.
%
%Examples
% > [total_error, max_error] = eq_area_error(2,10)
% total_error =
%    1.7764e-15
%  
% max_error =
%    4.4409e-16
%  
% > [total_error, max_error] = eq_area_error(3,1:6)
% total_error =
%    1.0e-12 *
%     0.0036    0.0036    0.1847    0.0142    0.0142    0.2132
%  
% max_error =
%    1.0e-12 *
%     0.0036    0.0018    0.1954    0.0284    0.0440    0.0777
%
%See also
% EQ_REGIONS, AREA_OF_SPHERE, AREA_OF_IDEAL_REGION

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
error(nargchk(2,2,nargin));
error(nargoutchk(2,2,nargout));

%
% Flatten N into a row vector.
%
shape = size(N);
n_partitions = prod(shape);
N = reshape(N,1,n_partitions);

total_error = zeros(size(N));
max_error   = zeros(size(N));
sphere_area = area_of_sphere(dim);

for partition_n = 1:n_partitions
    n = N(partition_n);
    regions = eq_regions(dim,n);
    ideal_area = area_of_ideal_region(dim,n);
    total_area = 0;
    for region_n = 1:size(regions,3)
        area = area_of_region(regions(:,:,region_n));
        total_area = total_area + area;
        region_error = abs(area - ideal_area);
        if region_error > max_error(partition_n)
            max_error(partition_n) = region_error;
        end
    end
    total_error(partition_n) = abs(sphere_area - total_area);
end
%
% Reshape output to same array size as original N.
%
total_error = reshape(total_error,shape);
max_error = reshape(max_error,shape);
%
% end function

function area = area_of_region(region)
%AREA_OF_REGION Area of given region
%
% area = area_of_region(region);

dim = size(region,1);
s_top = region(dim,1);
s_bot = region(dim,2);
if dim > 1
    area = area_of_collar(dim, s_top, s_bot)*area_of_region(region(1:dim-1,:))/area_of_sphere(dim-1);
else
    if s_bot == 0
        s_bot = 2*pi;
    end
    if s_top == s_bot
        s_bot = s_top + 2*pi;
    end
    area = s_bot - s_top;
end
%
% end function
