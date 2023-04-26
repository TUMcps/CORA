function dist = eq_min_dist(dim,N,varargin)
%EQ_MIN_DIST Minimum distance between center points of an EQ partition
%
%Syntax
% dist = eq_min_dist(dim,N,options);
%
%Description
% DIST = EQ_MIN_DIST(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region, and
% 3) sets DIST to be the minimum Euclidean distance between points of 
%    the EQ point set.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result DIST will be an array of the same size as N.
%
% DIST = EQ_MIN_DIST(dim,N,'offset','extra'), for dim == 2 or dim == 3, 
% uses exerimental extra rotation offsets to try to maximize the minimum 
% distance. For dim > 3, extra offsets are not used.
%
%Examples
% > dist=eq_min_dist(2,10)
%  dist =
%      1.0515
%  
% > dist=eq_min_dist(3,1:6)
%  dist =
%      2.0000    2.0000    1.4142    1.4142    1.4142    1.4142
%
%See also
% PARTITION_OPTIONS, EUCLIDEAN_DIST, EQ_ENERGY_DIST

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
dist = eq_point_set_property(@point_set_min_dist,dim,N,varargin{:});
%
% end function
