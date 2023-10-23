function res = test_zonoBundle_project
% test_zonoBundle_project - unit test function of project
%
% Syntax:
%    res = test_zonoBundle_project
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
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});
% compute projection (outer-approximation!)
projDim = 1;
zB_ = project(zB,projDim);

% sample points from original zonotope bundle
nrPoints = 10;
p = randPoint(zB,nrPoints,'standard');
% ensure that projected points lie inside computed projection
res(end+1,1) = all(contains(zB_,p(projDim,:)));

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
% compute projection (outer-approximation!)
projDim = 1;
zB_ = project(zB,projDim);

% sample points from original zonotope bundle
p = randPoint(zB,nrPoints,'standard');
% ensure that projected points lie inside computed projection
res(end+1,1) = all(contains(zB_,p(projDim,:)));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
