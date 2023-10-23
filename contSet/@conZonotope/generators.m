function G = generators(cZ)
% generators - returns the generators of a constrained zonotope as a matrix
%    whose column vectors are the generators
%
% Syntax:
%    G = generators(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    G - matrix of generators stored as column vectors
%
% Example: 
%    Z = [0 1 0 0 1;0 1 0 2 -1];
%    A = [-2 0 1 -1; 0 0 0 0]; b = [2;0];
%    cZ = conZonotope(Z,A,b);
%    G = generators(cZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       28-March-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

G = cZ.G;

% ------------------------------ END OF CODE ------------------------------
