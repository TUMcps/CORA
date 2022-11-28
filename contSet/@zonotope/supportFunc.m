function [val,x,fac] = supportFunc(Z,dir,varargin)
% supportFunc - calculates the upper or lower bound of a zonotope along a
%    certain direction
%
% Syntax:  
%    val = supportFunc(Z,dir)
%    val = supportFunc(Z,dir,type)
%
% Inputs:
%    Z - zonotope object
%    dir - direction for which the bounds are calculated (vector)
%    type - upper or lower bound ('lower' or 'upper')
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
% See also: conZonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_supportFunc('zonotope',Z,dir,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    val = vars{1}; x = []; fac = []; return
else
    Z = vars{1}; dir = vars{2}; type = vars{3};
end

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
    throw(CORAerror('CORAerror:notSupported',type));

end

% compute support vector
x = c + G*fac;

%------------- END OF CODE --------------