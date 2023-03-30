function [val,x,fac] = supportFunc_(Z,dir,type,varargin)
% supportFunc_ - calculates the upper or lower bound of a zonotope along a
%    certain direction
%
% Syntax:  
%    [val,x,fac] = supportFunc_(Z,dir)
%    [val,x,fac] = supportFunc_(Z,dir,type)
%
% Inputs:
%    Z - zonotope object
%    dir - direction for which the bounds are calculated (vector)
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the zonotope in the specified direction
%    x - support vector
%    fac - factor values that correspond to the upper bound
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/supportFunc_

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  10-December-2022 (MW, add type = 'range')
% Last revision:27-March-2023 (MW, rename supportFunc_)

%------------- BEGIN CODE --------------

% get object properties
c = center(Z);
G = generators(Z);

% project zonotope onto the direction
c_ = dir'*c;
G_ = dir'*G;

% upper or lower bound
if strcmp(type,'lower')
    val = c_ - sum(abs(G_));
    fac = -sign(G_)';
    
elseif strcmp(type,'upper')
    val = c_ + sum(abs(G_));
    fac = sign(G_)';
    
elseif strcmp(type,'range')
    val = interval(c_ - sum(abs(G_)), c_ + sum(abs(G_)));
    fac = [-sign(G_)' sign(G_)'];

end

% compute support vector
x = c + G*fac;

%------------- END OF CODE --------------