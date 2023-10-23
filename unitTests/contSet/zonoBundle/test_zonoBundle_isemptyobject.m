function res = test_zonoBundle_isemptyobject
% test_zonoBundle_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_zonoBundle_isemptyobject
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
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate zonotope bundles
zB1 = zonoBundle();

Z1 = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
Z2 = Z1 + [1;0;0];
zB2 = zonoBundle({Z1,Z2});

% check results
res = isemptyobject(zB1) && ~isemptyobject(zB2);

% ------------------------------ END OF CODE ------------------------------
