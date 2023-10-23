function res = test_emptySet_eq
% test_emptySet_eq - unit test function of eq
%
% Syntax:
%    res = test_emptySet_eq
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
res = O == O;

% different dimension
res(end+1,1) = ~(O == emptySet(n+1));

% init zonotope
Z = zonotope(zeros(2,1));
res(end+1,1) = ~(O == Z);

% init empty polytope
P = polytope([1 1;-1 -1],[1;-2]);
res(end+1,1) = eq(O,P);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
