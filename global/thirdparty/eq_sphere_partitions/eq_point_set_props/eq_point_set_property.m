function property = eq_point_set_property(fhandle,dim,N,varargin)
%EQ_POINT_SET_PROPERTY Property of an EQ point set
% 
%Syntax
% property = eq_point_set_property(fhandle,dim,N,options);
%
%Description
% PROPERTY = EQ_POINT_SET_PROPERTY(FHANDLE,dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) calls the function defined by FHANDLE which is expected to use the 
%    EQ point set to calculate the value of the result, PROPERTY.
%
% The argument FHANDLE must be a function handle. The function specified by
% FHANDLE must take as its argument a single (dim+1 by N) array, representing 
% N points of S^dim, in Cartesian coordinates, and must return a single value
% based on these points.
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result PROPERTY will be an array of the same size as N.
%
% PROPERTY = EQ_POINT_SET_PROPERTY(FHANDLE,dim,N,'offset','extra'), 
% for dim == 2 or dim == 3, uses exerimental extra rotation offsets to try to 
% maximize the minimum distance between points of the EQ point set.
% For dim > 3, extra offsets are not used.
%
%Examples
% > dist = eq_point_set_property(@point_set_min_dist,2,10)
%  dist =
%      1.0515
%
%See also
% EQ_POINT_SET, FEVAL, PARTITION_OPTIONS

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
error(nargchk(3,5,nargin));
%
% dim is the number of dimensions
% N is the number of regions
%
% If the option 'offset' is 'extra', then use experimental extra offsets
% for S^2 and S^3 regions to try to minimize energy.
% The default is not to use extra offsets.
%
pdefault.extra_offset =  false;
popt = partition_options(pdefault, varargin{:});
extra_offset = popt.extra_offset;
%
% Flatten N into a row vector.
%
shape = size(N);
n_partitions = prod(shape);
N = reshape(N,1,n_partitions);

property = zeros(size(N));
for partition_n = 1:n_partitions
    points = eq_point_set(dim,N(partition_n),extra_offset);
    property(partition_n) = feval(fhandle,points);
end    
%
% Reshape output to same array size as original N.
%
property = reshape(property,shape);
%
% end function
