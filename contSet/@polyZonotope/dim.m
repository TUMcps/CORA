function n = dim(pZ)
% dim - returns the dimension of the ambient space of a polynomial zonotope
%
% Syntax:
%    n = dim(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    GI = [0;0.5];
%    E = [1 0 3;0 1 1];
%    pZ = polyZonotope(c,G,GI,E);
%    n = dim(pZ)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m

% Authors:       Niklas Kochdumper
% Written:       26-December-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = size(pZ.c,1);

% ------------------------------ END OF CODE ------------------------------
