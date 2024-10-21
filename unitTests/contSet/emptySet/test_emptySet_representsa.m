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

% init empty set
n = 2;
O = emptySet(n);

% compare to other representations
assert(~representsa(O,'origin'));
assert(~representsa(O,'point'));

assert(representsa(O,'emptySet'));
assert(representsa(O,'interval'));
assert(representsa(O,'zonotope'));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
