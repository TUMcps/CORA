function cZ = conZonotope(Z)
% conZonotope - convert a zonotope object into a conZonotope object
%
% Syntax:
%    cZ = conZonotope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    cZ - conZonotope object
%
% Example:
%    Z = zonotope([1;0],[1 0; 0 1]);
%    cZ = conZonotope(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       23-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call constructor with center and generator matrix as inputs
cZ = conZonotope(Z.c, Z.G);

% ------------------------------ END OF CODE ------------------------------
