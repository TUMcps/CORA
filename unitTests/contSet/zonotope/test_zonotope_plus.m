function res = test_zonotope_plus
% test_zonotope_plus - unit test function of plus
%
% Syntax:
%    res = test_zonotope_plus
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
% See also: -

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D zonotopes
Z1 = zonotope([-4; 1], [-3, -2, -1; 2, 3, 4]);
Z2 = zonotope([1;-1], [10; -10]);
Z_plus = Z1 + Z2;
% compare to true result
c_plus = [-3; 0];
G_plus = [-3, -2, -1, 10; ...
            2, 3, 4, -10];
assert(compareMatrices(Z_plus.c,c_plus));
assert(compareMatrices(Z_plus.G,G_plus));

% Minkowski sum with empty set
Z_empty = zonotope.empty(2);
Z_plus = Z1 + Z_empty;
assert(representsa(Z_plus,'emptySet'));


% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
