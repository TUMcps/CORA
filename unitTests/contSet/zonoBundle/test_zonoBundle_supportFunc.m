function res = test_zonoBundle_supportFunc
% test_zonoBundle_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = test_zonoBundle_supportFunc
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

res = true;
% fully-empty set
zB = zonoBundle.empty(2);
[val,x] = supportFunc(zB,[1;0],'upper');
res(end+1,1) = val == -Inf && isnumeric(x) && isempty(x);
[val,x] = supportFunc(zB,[1;0],'lower');
res(end+1,1) = val == Inf && isnumeric(x) && isempty(x);
val = supportFunc(zB,[1;0],'range');
res(end+1,1) = isequal(val,interval(-Inf,Inf));

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});

dir = [1;0];
[val,x] = supportFunc(zB,dir,'upper');
res(end+1,1) = withinTol(val,4) && compareMatrices(x,[4;0]);
val = supportFunc(zB,dir,'lower');
res(end+1,1) = withinTol(val,-2);
dir = [-1;-1];
[~,x] = supportFunc(zB,dir,'upper');
res(end+1,1) = compareMatrices(x,[-2;-1]);
dir = [0,1];
val = supportFunc(zB,dir,'lower');
res(end+1,1) = withinTol(val,-1);

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
val = supportFunc(zB,[1;0]);
res(end+1,1) = val == -Inf;
val = supportFunc(zB,[1;0],'lower');
res(end+1,1) = val == Inf;
val = supportFunc(zB,[1;0],'range');
res(end+1,1) = isequal(val,interval(-Inf,Inf));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
