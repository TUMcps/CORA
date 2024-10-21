function res = test_polyZonotope_mtimes
% test_polyZonotope_mtimes - unit test function for multiplication between
%    a matrix and a polynomial zonotope
%
% Syntax:
%    res = test_polyZonotope_mtimes
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       26-June-2018
% Last update:   05-October-2024 (MW, test projections, regular matrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D, matrix * empty zonotope
pZ = polyZonotope.empty(2);
% ...square matrix
M = [-1 2; 3 -4];
pZ_mtimes = M*pZ;
assert(representsa(pZ_mtimes,'emptySet') && dim(pZ_mtimes) == 2);
% ...projection onto subspace
M = [-1 2];
pZ_mtimes = M*pZ;
assert(representsa(pZ_mtimes,'emptySet') && dim(pZ_mtimes) == 1);
% ...projection to higher-dimensional space
M = [-1 2; 3 -4; 5 6];
pZ_mtimes = M*pZ;
assert(representsa(pZ_mtimes,'emptySet') && dim(pZ_mtimes) == 3);

% 2D, create polynomial zonotope
c = [1; 2]; G = [1 2 1 -3; 1 -1 2 -1]; E = [1 0 0 2; 0 1 2 1]; GI = [0 -1; 1 1];
pZ = polyZonotope(c,G,GI,E);
% ...square matrix
M = [1 -0.5; -1 0];
pZ_mtimes = M * pZ;
assert(compareMatrices(pZ_mtimes.c,M*c));
assert(compareMatrices(pZ_mtimes.G,M*G));
assert(compareMatrices(pZ_mtimes.GI,M*GI));
assert(compareMatrices(pZ_mtimes.E,pZ.E));
assert(compareMatrices(pZ_mtimes.id,pZ.id));
% ...projection onto subspace
M = [-1 2];
pZ_mtimes = M*pZ;
assert(compareMatrices(pZ_mtimes.c,M*c));
assert(compareMatrices(pZ_mtimes.G,M*G));
assert(compareMatrices(pZ_mtimes.GI,M*GI));
assert(compareMatrices(pZ_mtimes.E,pZ.E));
assert(compareMatrices(pZ_mtimes.id,pZ.id));
% ...projection to higher-dimensional space
M = [-1 2; 3 -4; 5 6];
pZ_mtimes = M*pZ;
assert(compareMatrices(pZ_mtimes.c,M*c));
assert(compareMatrices(pZ_mtimes.G,M*G));
assert(compareMatrices(pZ_mtimes.GI,M*GI));
assert(compareMatrices(pZ_mtimes.E,pZ.E));
assert(compareMatrices(pZ_mtimes.id,pZ.id));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
