function res = test_emptySet_ne
% test_emptySet_ne - unit test function of ne
%
% Syntax:
%    res = test_emptySet_ne
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
assert(~(O ~= O));

% different dimension
assert(O ~= emptySet(n+1));

% init zonotope
Z = zonotope(zeros(2,1));
assert(O ~= Z);

% init empty polytope
P = polytope([1 1;-1 -1],[1;-2]);
assert(~ne(O,P));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
