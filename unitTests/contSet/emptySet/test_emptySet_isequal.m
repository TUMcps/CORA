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
res = isequal(O,O);

% different dimension
res(end+1,1) = ~isequal(O,emptySet(n+1));

% init zonotope
Z = zonotope(zeros(2,1));
res(end+1,1) = ~isequal(O,Z);

% init empty polytope
P = polytope([1 1;-1 -1],[1;-2]);
res(end+1,1) = isequal(O,P);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
