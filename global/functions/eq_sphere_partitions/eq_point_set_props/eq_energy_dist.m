function [energy,dist] = eq_energy_dist(dim,N,s,varargin)
%EQ_ENERGY_DIST Energy and minimum distance of an EQ point set
%
% Syntax
% [energy,dist] = eq_energy_dist(dim,N,s,options);
%
%Description
% [ENERGY,DIST] = EQ_ENERGY_DIST(dim,N,s) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to 
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) sets ENERGY to be the r^(-s) energy of the EQ point set, and
% 4) optionally, sets DIST to be the minimum Euclidean distance between 
%    points of the EQ point set.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The results ENERGY and DIST will be arrays of the same size as N.
% The result DIST is optional.
%
% [ENERGY,DIST] = EQ_ENERGY_DIST(dim,N) uses the default value dim-1 for s.
%
% [ENERGY,DIST] = EQ_ENERGY_DIST(dim,N,s,'offset','extra') uses experimental
% extra offsets for S^2 and S^3 to try to minimize energy. 
% For dim > 3, extra offsets are not used.
%
%Examples
% > energy=eq_energy_dist(2,10)
%  energy =
%     32.7312
%
% > [energy,dist]=eq_energy_dist(3,1:6,0)
%  energy =
%           0   -0.6931   -1.3863   -2.7726   -4.1589   -6.2383
%  dist =
%      2.0000    2.0000    1.4142    1.4142    1.4142    1.4142
%  
% > [energy,dist]=eq_energy_dist(3,100,1,'offset','extra')
%  energy =
%     4.0042e+03
%  dist =
%      0.6545
%
%See also
% EQ_POINT_SET, PARTITION_OPTIONS, POINT_SET_ENERGY_DIST, EUCLIDEAN_DIST, 
% MIN_DIST

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments.
%
error(nargchk(2,5,nargin));
%
if nargin < 3
    %
    % The default value of s is dim-1.
    %
    s = dim-1;
elseif ischar(s)
    %
    % Ensure that the user has not omitted argument s.
    %
    error('Argument s must be numeric.');    
end
%
% If the option 'offset' is 'extra', then use experimental extra offsets
% for S^2 and S^3 regions to try to minimize energy.
% The default is not to use extra offsets.
%
pdefault.extra_offset =  false;
if nargin < 4
    extra_offset = pdefault.extra_offset;
else
    popt = partition_options(pdefault, varargin{:});
    extra_offset = popt.extra_offset;
end
%
% Flatten N into a row vector.
%
shape = size(N);
n_partitions = prod(shape);
N = reshape(N,1,n_partitions);

energy = zeros(size(N));
if nargout > 1
    dist = zeros(size(N));
end
for partition_n = 1:n_partitions
    points = eq_point_set(dim,N(partition_n),extra_offset);
    %
    if nargout > 1
        [energy(partition_n),dist(partition_n)] = point_set_energy_dist(points,s);
    else
        energy(partition_n) = point_set_energy_dist(points,s);
    end
end
%
% Reshape output to same array size as original N.
%
energy = reshape(energy,shape);
if nargout > 1
    dist = reshape(dist,shape);
end    
%
% end function
