function res = test_zonoBundle_copy
% test_zonoBundle_copy - unit test function of copy
%
% Syntax:
%    res = test_zonoBundle_copy
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
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D zonoBundle
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
zB_copy = copy(zB);
assert(isequal(zB,zB_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
