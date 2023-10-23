function res = isFullDim(cPZ)
% isFullDim - checks if the dimension of the affine hull of a constrained
%    polynomial zonotope is equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    res - true/false
%
% Example:
%    cPZ1 = conPolyZono([1;3],[1 2;1 -2],[1 2;1 1],[1;0]);
%    cPZ2 = conPolyZono([1;3],[1 2;1 2],[1 2;1 1],[-1;-1]);
%
%    isFullDim(cPZ1)
%    isFullDim(cPZ2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/isFullDim

% Authors:       Niklas Kochdumper
% Written:       08-February-2020 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get dimension and rank
if representsa_(cPZ,'emptySet',eps)
    res = false; return;
end
n = length(cPZ.c);
rk = rank([cPZ.G,cPZ.GI]);

% compare dimension and rank
res = n == rk;

% ------------------------------ END OF CODE ------------------------------
