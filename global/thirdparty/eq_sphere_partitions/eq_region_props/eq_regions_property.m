function property = eq_regions_property(fhandle,dim,N)
%EQ_REGIONS_PROPERTY Property of regions of an EQ partition
% 
%Syntax
% property = eq_regions_property(fhandle,dim,N);
%
%Description
% PROPERTY = EQ_REGIONS_PROPERTY(FHANDLE,dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) calls the function defined by FHANDLE which is expected to use the 
%    regions to calculate the value of the result, PROPERTY.
%
% The argument FHANDLE must be a function handle. The function specified by
% FHANDLE must take as its argument a single (dim by 2 by N) array, 
% representing N regions of S^dim, in spherical polar coordinates, and must 
% return a single value based on these regions.
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result PROPERTY will be an array of the same size as N.
%
%Examples
% See code in Matlab M files eq_diam_bound.m, eq_vertex_diam.m.
%
%See also 
% EQ_REGIONS, FEVAL, EQ_DIAM_BOUND, EQ_VERTEX_DIAM

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Flatten N into a row vector.
%
shape = size(N);
n_partitions = prod(shape);
N = reshape(N,1,n_partitions);

property = zeros(size(N));
for partition_n = 1:n_partitions
    regions = eq_regions(dim,N(partition_n));
    property(partition_n) = feval(fhandle,regions);
end    
%
% Reshape output to same array size as original N.
%
property = reshape(property,shape);
%
% end function
