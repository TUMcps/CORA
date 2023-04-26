function res = test_zonoBundle_interval
% test_zonoBundle_interval - unit test function of interval conversion
%
% Syntax:  
%    res = test_zonoBundle_interval
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
I = interval(zB);
res = isempty(I);

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
% convert to interval
I = interval(zB);
res(end+1,1) = contains(I,zB);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
% convert to interval
I = interval(zB);
res(end+1,1) = isempty(I);

% combine results
res = all(res);

%------------- END OF CODE --------------
