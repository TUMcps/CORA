function res = test_emptySet_representsa
% test_emptySet_representsa - unit test function of representsa
%
% Syntax:
%    res = test_emptySet_representsa
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
% Written:       25-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% init empty set
n = 2;
O = emptySet(n);

% compare to other representations
res(end+1,1) = ~representsa(O,'origin');
res(end+1,1) = ~representsa(O,'point');

res(end+1,1) = representsa(O,'emptySet');
res(end+1,1) = representsa(O,'interval');
res(end+1,1) = representsa(O,'zonotope');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
