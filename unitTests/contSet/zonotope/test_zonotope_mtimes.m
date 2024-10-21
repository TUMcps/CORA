function res = test_zonotope_mtimes
% test_zonotope_mtimes - unit test function of mtimes
%
% Syntax:
%    res = test_zonotope_mtimes
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   05-October-2024 (MW, test projections)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D, matrix * empty zonotope
Z = zonotope.empty(2);
% ...square matrix
M = [-1 2; 3 -4];
Z_mtimes = M*Z;
assert(representsa(Z_mtimes,'emptySet') && dim(Z_mtimes) == 2);
% ...projection onto subspace
M = [-1 2];
Z_mtimes = M*Z;
assert(representsa(Z_mtimes,'emptySet') && dim(Z_mtimes) == 1);
% ...projection to higher-dimensional space
M = [-1 2; 3 -4; 5 6];
Z_mtimes = M*Z;
assert(representsa(Z_mtimes,'emptySet') && dim(Z_mtimes) == 3);

% 2D, matrix * non-empty zonotope
Z = zonotope([-4;1], [-3 -2 -1; 2 3 4]);
% ...square matrix
M = [-1 2; 3 -4];
Z_mtimes = M*Z;
Z_true = zonotope([6; -16], [7, 8, 9; -17, -18, -19]);
assert(isequal(Z_mtimes,Z_true));
% ...projection onto subspace
M = [-1 2];
Z_mtimes = M*Z;
Z_true = zonotope(6,24);
assert(isequal(Z_mtimes,Z_true));
% ...projection to higher-dimensional space
M = [-1 2; 1 0; 0 1];
Z_mtimes = M*Z;
Z_true = zonotope([6;-4;1],[7 8 9; -3 -2 -1; 2 3 4]);
assert(isequal(Z_mtimes,Z_true));

% 2D, non-empty zonotope * matrix/scalar
Z = zonotope([-4;1], [-3 -2 -1; 2 3 4]);
M = 2;
Z_mtimes = Z * M;
Z_true = zonotope([-8;2], [-6 -4 -2; 4 6 8]);
assert(isequal(Z_mtimes,Z_true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
