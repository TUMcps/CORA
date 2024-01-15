function res = test_zonoBundle_representsa
% test_zonoBundle_representsa - unit test function of representsa
%
% Syntax:
%    res = test_zonoBundle_representsa
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: 20-July-2023 (MW, rename '...representsa')

% ------------------------------ BEGIN CODE -------------------------------

% fully-empty zonoBundle
zB = zonoBundle.empty(2);
res = representsa(zB,'emptySet');

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
res(end+1,1) = ~representsa(zB,'emptySet');

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
res(end+1,1) = representsa(zB,'emptySet');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
