function [bound_coeff,vertex_coeff] =  eq_diam_coeff(dim,N)
%EQ_DIAM_COEFF Coefficients of diameter bound and vertex diameter of EQ partition
%
%Syntax
% [bound_coeff,vertex_coeff] = eq_diam_coeff(dim,N);
%
%Description
% [BOUND_COEFF,VERTEX_COEFF] = EQ_DIAM_COEFF(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) finds the maximum of the per-region diameter bound over all the regions 
%    of the partition,
% 3) sets BOUND_COEFF to be the diameter bound coefficient, defined as the
%    solution to
% 
%    max_diam_bound == BOUND_COEFF N^(-1/dim),
%
% 4) optionally finds the maximum vertex diameter over all the regions of the
%    partition, and
% 5) optionally sets VERTEX_COEFF to be the vertex diameter coefficient,
%    defined as the solution to
% 
%    max_vertex_diam == VERTEX_COEFF N^(-1/dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result BOUND_COEFF and the optional result VERTEX_COEFF will be arrays of
% the same size as N.
%
%Examples
% > bound_coeff=eq_diam_coeff(2,10)
%  bound_coeff =
%      5.2915
%  
% > [bound_coeff,vertex_coeff]=eq_diam_coeff(3,1:6)
%  bound_coeff =
%      2.0000    2.5198    2.8845    3.1748    3.4200    3.6342
%  
%  vertex_coeff =
%      2.0000    2.5198    2.8845    3.1748    3.4200    3.6342
%
%See also 
% EQ_DIAM_BOUND, EQ_VERTEX_DIAM, EQ_REGIONS, EQ_VERTEX_DIAM_COEFF
 
% Copyright 2004-2005 Paul Leopardi for the University of NSW.
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
error(nargoutchk(0,2,nargout));

if nargout < 2
    bound_coeff = eq_diam_bound(dim,N) .* N.^(1/dim);
else
    %
    % Flatten N into a row vector.
    %
    shape = size(N);
    n_partitions = prod(shape);
    N = reshape(N,1,n_partitions);
    
    bound_coeff =  zeros(size(N));
    vertex_coeff = zeros(size(N));
    for partition_n = 1:n_partitions
        n = N(partition_n);
        regions = eq_regions(dim,n);
        scale = n^(1/dim);
        bound_coeff(partition_n) =  max_diam_bound_of_regions(regions)  * scale;
        vertex_coeff(partition_n) = max_vertex_diam_of_regions(regions) * scale;
    end    
    %
    % Reshape output to same array size as original N.
    %
    bound_coeff =  reshape(bound_coeff,shape);
    vertex_coeff = reshape(vertex_coeff,shape);
end
%
%end function
