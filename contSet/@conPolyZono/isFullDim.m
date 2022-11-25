function res = isFullDim(cPZ)
% isFullDim - check if a constrained polynomial zonotope is 
%             full-dimensional
%
% Syntax:  
%    res = isFullDim(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    res - 1 if conPolyZono is full-dimensional, 0 else
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

% Author:       Niklas Kochdumper
% Written:      08-February-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % get dimension and rank
    n = length(cPZ.c);
    rk = rank([cPZ.G,cPZ.Grest]);

    % compare dimension and rank
    res = n == rk;

%------------- END OF CODE --------------