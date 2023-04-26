function res = test_zonoBundle_isempty
% test_zonoBundle_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_zonoBundle_isempty
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

% Author:       Mark Wetzlinger
% Written:      23-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% fully-empty zonoBundle
zB = zonoBundle();
res = isempty(zB);

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
res(end+1,1) = ~isempty(zB);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
res(end+1,1) = isempty(zB);

% combine results
res = all(res);

%------------- END OF CODE --------------
