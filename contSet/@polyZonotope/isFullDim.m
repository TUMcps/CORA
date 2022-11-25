function res = isFullDim(pZ)
% isFullDim - check if a polynomial zonotope is full-dimensional
%
% Syntax:  
%    res = isFullDim(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - 1 if polyZonotope is full-dimensional, 0 else
%
% Example:
%    pZ1 = polyZonotope([1;3],[1 2;1 -2],[1;0],[1 2;1 1]);
%    pZ2 = polyZonotope([1;3],[1 2;1 2],[-1;-1],[1 2;1 1]);
%
%    isFullDim(pZ1)
%    isFullDim(pZ2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Author:       Niklas Kochdumper
% Written:      02-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% get dimension and rank
n = length(pZ.c);
rk = rank([pZ.G,pZ.Grest]);

% compare dimension and rank
res = n == rk;

%------------- END OF CODE --------------