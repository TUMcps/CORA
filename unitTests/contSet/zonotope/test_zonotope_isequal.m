function res = test_zonotope_isequal
% test_zonotope_isequal - unit test function of isequal
%
% Syntax:
%    res = test_zonotope_isequal
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

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   21-April-2020
%                09-August-2020 (enhance randomness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D, minimal and non-minimal
Z1 = zonotope(4,[1 3 -2]);
Z2 = zonotope(4,6);
assert(isequal(Z1,Z2));
assert(isequal(Z2,Z1));


% 2D, different order of generators
Z1 = zonotope(ones(2,1),[1 2 5 3 3;
                          2 3 0 4 1]);
Z2 = zonotope(ones(2,1),[2 1 3 5 3;
                          3 2 4 0 1]);
assert(isequal(Z1,Z2));
assert(isequal(Z2,Z1));


% 2D different sign
Z1 = zonotope([0;1],[1 0 -1; 1 1 2]);
Z2 = zonotope([0;1],[-1 0 -1; -1 -1 2]);
assert(isequal(Z1,Z2));
assert(isequal(Z2,Z1));


% 3D, with zero-generator
Z1 = zonotope([1;5;-1],[2 4; 6 0; 4 8]);
Z2 = zonotope([1;4;-1],[2 4; 6 0; 4 8]);
Z3 = zonotope([1;5;-1],[2 0 4; 6 0 0; 4 0 8]);
assert(isequal(Z1,Z3));
assert(isequal(Z3,Z1));
assert(~isequal(Z1,Z2));
assert(~isequal(Z2,Z1));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
