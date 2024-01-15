function res = test_zonoBundle_isFullDim
% test_zonoBundle_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_zonoBundle_isFullDim
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% fully-empty zonoBundle
zB = zonoBundle.empty(2);
res = ~isFullDim(zB);

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
res(end+1,1) = isFullDim(zB);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
res(end+1,1) = ~isFullDim(zB);

% intersection of 2D sets is just a line
Z1 = zonotope([2;0],[1 0; 0 1]);
Z2 = zonotope([0;0],[1 0; 0 1]);
zB = zonoBundle({Z1,Z2});
res(end+1,1) = ~isFullDim(zB);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
