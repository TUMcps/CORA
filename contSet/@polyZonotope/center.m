function c = center(pZ)
% center - returns the starting point of a polynomial zonotope, which does
%    not correspond to the geometric center of the polynomial zonotope
%    in general
%
% Syntax:
%    c = center(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    c - starting point of the polynomial zonotope pZ
%
% Example: 
%    pZ = polyZonotope([3;4],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    center(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/center, polyZonotope/center

% Authors:       Niklas Kochdumper
% Written:       27-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = pZ.c;

% ------------------------------ END OF CODE ------------------------------
