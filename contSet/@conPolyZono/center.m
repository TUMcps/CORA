function c = center(cPZ)
% center - returns the starting point of a constrained polynomial zonotope
%
% Syntax:
%    c = center(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    c - starting point of the constrained polynomial zonotope
%
% Example:
%    A = 1/8 * [-10 2 2 3 3];
%    b = -3/8;
%    EC = [1 0 1 2 0; 0 1 1 0 2; 0 0 0 0 0];
%    c = [0;0];
%    G = [1 0 1 -1/4;0 1 -1 1/4];
%    E = [1 0 2 0;0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%    c = center(cPZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/center, polyZonotope/center

% Authors:       Niklas Kochdumper
% Written:       14-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
                     
c = cPZ.c;

% ------------------------------ END OF CODE ------------------------------
