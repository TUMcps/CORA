function vertex_diam = eq_vertex_diam(dim,N)
%EQ_VERTEX_DIAM Maximum vertex diameter of EQ partition
%
%Syntax
% vertex_diam = eq_vertex_diam(dim,N);
%
%Description
% VERTEX_DIAM = EQ_VERTEX_DIAM(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) sets VERTEX_DIAM to be the maximum vertex diameter over all the regions
%    of the partition.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result VERTEX_DIAM will an array of the same size as N.
%
%Examples
% > vertex_diam=eq_vertex_diam(2,10)
%  vertex_diam =
%      1.4142
%  
% > vertex_diam=eq_vertex_diam(3,1:6)
%  vertex_diam =
%       2     2     2     2     2     2
%  
%See also
% EQ_DIAM_BOUND, EQ_VERTEX_DIAM_COEFF, EQ_DIAM_COEFF, EQ_REGIONS_PROPERTY

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

vertex_diam = eq_regions_property(@max_vertex_diam_of_regions,dim,N);
%
% end function
