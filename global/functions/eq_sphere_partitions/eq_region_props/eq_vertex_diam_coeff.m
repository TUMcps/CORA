function coeff =  eq_vertex_diam_coeff(dim,N)
%EQ_VERTEX_DIAM_COEFF Coefficient of maximum vertex diameter of EQ partition
%
%Syntax
% coeff = eq_vertex_diam_coeff(dim,N);
%
%Description
% COEFF = EQ_VERTEX_DIAM_COEFF(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) finds the maximum vertex diameter over all the regions of the partition,
% 3) sets COEFF to be the vertex diameter coefficient,
%    defined as the solution to
% 
%    max_vertex_diam == COEFF N^(-1/dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result COEFF will an array of the same size as N.
%
%Examples
% > coeff=eq_vertex_diam_coeff(2,10)
%  coeff =
%      4.4721
%  
% > coeff=eq_vertex_diam_coeff(3,1:6)
%  coeff =
%      2.0000    2.5198    2.8845    3.1748    3.4200    3.6342
%
%See also
% EQ_VERTEX_DIAM, EQ_DIAM_COEFF

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

coeff = eq_vertex_diam(dim,N) .* N.^(1/dim);
%
%end function