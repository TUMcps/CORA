function res = test_zonoBundle_randPoint
% test_zonoBundle_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_zonoBundle_randPoint
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
p = randPoint(zB);
assert(isnumeric(p) && isempty(p) && size(p,1) == n);

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});

% sample random points with different syntax
p = randPoint(zB);
p = [p, randPoint(zB,5)];
p = [p, randPoint(zB,5,'extreme')];
p = [p, randPoint(zB,'all')];

% check if all points are contained
assert(all(contains(zB,p,'exact',1e-12)));

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
% sample random point
p = randPoint(zB);
assert(isnumeric(p) && isempty(p));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
