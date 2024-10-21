function res = test_emptySet_isequal
% test_emptySet_isequal - unit test function of isequal
%
% Syntax:
%    res = test_emptySet_isequal
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

% Authors:       Mark Wetzlinger
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init empty set
n = 2;
O = emptySet(n);

% comparison with itself
assert(isequal(O,O));

% different dimension
assert(~isequal(O,emptySet(n+1)));

% init zonotope
Z = zonotope(zeros(2,1));
assert(~isequal(O,Z));

% init empty polytope
P = polytope([1 1;-1 -1],[1;-2]);
assert(isequal(O,P));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
