function res = test_zonoBundle_center
% test_zonoBundle_center - unit test function of center
%
% Syntax:
%    res = test_zonoBundle_center
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
n = 2;
zB = zonoBundle.empty(n);
c = center(zB);
res = isnumeric(c) && isempty(c) && size(c,1) == n;

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
% compute center
c = center(zB);
res(end+1,1) = contains(zB,c);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
% compute center
c = center(zB);
res(end+1,1) = isnumeric(c) && isempty(c);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
