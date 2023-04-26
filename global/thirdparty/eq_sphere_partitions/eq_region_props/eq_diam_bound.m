function diam_bound = eq_diam_bound(dim,N)
%EQ_DIAM_BOUND Maximum per-region diameter bound of EQ partition
%
%Syntax
% diam_bound = eq_diam_bound(dim,N);
%
%Description
% DIAM_BOUND = EQ_DIAM_BOUND(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) sets DIAM_BOUND to be the maximum of the per-region diameter bound over
%    all the regions of the partition.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result DIAM_BOUND will an array of the same size as N.
%
%Examples
% > diam_bound=eq_diam_bound(2,10)
%  diam_bound =
%      1.6733
%  
% > diam_bound=eq_diam_bound(3,1:6)
%  diam_bound =
%       2     2     2     2     2     2
%  
%See also
% EQ_VERTEX_DIAM, EQ_DIAM_COEFF, EQ_REGIONS_PROPERTY

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

diam_bound = eq_regions_property(@max_diam_bound_of_regions,dim,N);
%
% end function
