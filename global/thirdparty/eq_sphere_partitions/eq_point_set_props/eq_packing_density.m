function density = eq_packing_density(dim,N,varargin)
%EQ_PACKING_DENSITY  Density of packing given by minimum distance of EQ point set
%
%Syntax
% density = eq_packing_density(dim,N,options);
%
%Description
% DENSITY = EQ_PACKING_DENSITY(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) finds the minimum Euclidean distance between points of the EQ point set,
% 4) sets DENSITY to be the maximum density of a packing of S^dim by equal 
%    spherical caps with centers at the EQ point set.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result DENSITY will be an array of the same size as N.
%
% DENSITY = EQ_PACKING_DENSITY(dim,N,'offset','extra'), for dim == 2 or dim == 3, 
% uses experimental extra rotation offsets to try to maximize the minimum 
% distance. For dim > 3, extra offsets are not used.
%
%Notes
% The packing density is defined to be the sum of the areas of the spherical 
% caps of the packing, divided by the area of the unit sphere S^dim.
%
% The spherical radius of the caps in the packing is half the minimum spherical
% distance between center points. The spherical radius for N == 1 is a special 
% case. It is defined to be pi.
%
%Examples 
% > density=eq_packing_density(2,10)
%  density =
%      0.7467
%  
% > density=eq_packing_density(3,1:6)
%  density =
%      1.0000    1.0000    0.2725    0.3634    0.4542    0.5451
%
%See also
% EQ_MIN_DIST, AREA_OF_CAP, AREA_OF_SPHERE, PARTITION_OPTIONS

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
error(nargchk(2,4,nargin));
%
min_euclidean_dist = eq_min_dist(dim,N,varargin{:});
density = calc_packing_density(dim,N,min_euclidean_dist);
%
% end function
