function res = isFullDim(Z)
% isFullDim - check if a zonotope is full-dimensional
%
% Syntax:  
%    res = isFullDim(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    res - 1 if zonotope is full-dimensional, 0 else
%
% Example:
%    zono1 = zonotope([1 2 1;3 1 2]);
%    zono2 = zonotope([1 2 1;3 4 2]);
%
%    isFullDim(zono1)
%    isFullDim(zono2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isempty

% Author:       Niklas Kochdumper
% Written:      02-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    Zdim = dim(Z);
    Grank = rank(generators(Z));
    
    res = Zdim == Grank;

%------------- END OF CODE --------------