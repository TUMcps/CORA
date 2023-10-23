function res = isFullDim(pZ)
% isFullDim - checks if the dimension of the affine hull of a polynomial
%    zonotope is equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - true/false
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

% Authors:       Niklas Kochdumper
% Written:       02-January-2020 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if representsa_(pZ,'emptySet',eps)
    res = false; return;
end

% get dimension and rank
n = length(pZ.c);
rk = rank([pZ.G,pZ.GI]);

% compare dimension and rank
res = n == rk;

% ------------------------------ END OF CODE ------------------------------
