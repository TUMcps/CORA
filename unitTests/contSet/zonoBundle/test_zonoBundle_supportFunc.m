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

tol = 1e-10;

% fully-empty set
zB = zonoBundle.empty(2);
[val,x] = supportFunc(zB,[1;0],'upper');
assert(val == -Inf && isnumeric(x) && isempty(x));
[val,x] = supportFunc(zB,[1;0],'lower');
assert(val == Inf && isnumeric(x) && isempty(x));
val = supportFunc(zB,[1;0],'range');
assert(isequal(val,interval(-Inf,Inf)));

% non-empty intersection
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});

dir = [1;0];
[val,x] = supportFunc(zB,dir,'upper');
assert(withinTol(val,4) && compareMatrices(x,[4;0],tol));
val = supportFunc(zB,dir,'lower');
assert(withinTol(val,-2));
dir = [-1;-1];
[~,x] = supportFunc(zB,dir,'upper');
assert(compareMatrices(x,[-2;-1],tol));
dir = [0,1];
val = supportFunc(zB,dir,'lower');
assert(withinTol(val,-1));

% empty intersection
Z2 = zonotope([-4;1],[0.5 1; 1 -1]);
zB = zonoBundle({Z1,Z2});
val = supportFunc(zB,[1;0]);
assert(val == -Inf);
val = supportFunc(zB,[1;0],'lower');
assert(val == Inf);
val = supportFunc(zB,[1;0],'range');
assert(isequal(val,interval(-Inf,Inf)));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
