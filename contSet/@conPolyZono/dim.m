function n = dim(cPZ)
% dim - returns the dimension of the ambient space of the constrained
%    polynomial zonotope
%
% Syntax:
%    n = dim(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    A = 1/8 * [-10 2 2 3 3];
%    b = -3/8;
%    EC = [1 0 1 2 0; 0 1 1 0 2; 0 0 0 0 0];
%    c = [0;0];
%    G = [1 0 1 -1/4;0 1 -1 1/4];
%    E = [1 0 2 0;0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%    n = dim(cPZ)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/dim

% Authors:       Niklas Kochdumper
% Written:       21-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = size(cPZ.c,1);

% ------------------------------ END OF CODE ------------------------------
