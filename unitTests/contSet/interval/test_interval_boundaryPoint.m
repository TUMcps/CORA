function res = test_interval_boundaryPoint
% test_interval_boundaryPoint - unit test function of boundary point
%    computation for intervals
%
% Syntax:
%    res = test_interval_boundaryPoint
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

% Authors:       Mark Wetzlinger
% Written:       13-August-2024
% Last update:   17-October-2024 (MW, test for start point)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-14;

% empty set
I = interval.empty(2);
dir = [1;1];
x = boundaryPoint(I,dir);
assert(all(size(x) == [2,0]));

% bounded, non-degenerate
I = interval([-2;1],[4;2]);
dir = [1;1];
x = boundaryPoint(I,dir);
x_true = [1.5;2];
assert(compareMatrices(x,x_true,tol,'equal',true));

% bounded, non-degenerate, different start point
I = interval([-2;1],[4;2]);
dir = [1;1];
startPoint = [-1;1];
x = boundaryPoint(I,dir,startPoint);
x_true = [0;2];
assert(compareMatrices(x,x_true,tol,'equal',true));

% bounded, degenerate
I = interval([-1;0],[-1;4]);
dir = [1;1];
x = boundaryPoint(I,dir);
x_true = center(I);
assert(compareMatrices(x,x_true,tol,'equal',true));

% bounded, degenerate, different start point
I = interval([-1;0],[-1;4]);
dir = [1;1];
startPoint = [-1;0];
x = boundaryPoint(I,dir,startPoint);
assert(compareMatrices(x,startPoint,tol,'equal',true));

% unbounded
I = interval([-Inf; 1],[0;2]);
dir = [-1;1];
startPoint = [0;1];
x = boundaryPoint(I,dir,startPoint);
assert(all(x == [-Inf;2]));

% unbounded, but vector does reach a finite boundary
I = interval([-Inf; 1],[0;2]);
dir = [0;1];
startPoint = [0;1];
x = boundaryPoint(I,dir,startPoint);
x_true = [0;2];
assert(compareMatrices(x,x_true,tol,'equal',true));


% wrong calls
I = interval([-2;1],[4;2]);
% ...all-zero 'direction'
assertThrowsAs(@boundaryPoint,'CORA:wrongValue',I,[0;0]);
% ...start point not in the set
assertThrowsAs(@boundaryPoint,'CORA:wrongValue',I,[1;1],[-5;10]);
% ...dimension mismatch
assertThrowsAs(@boundaryPoint,'CORA:dimensionMismatch',I,[1;1;1],[-5;10]);
assertThrowsAs(@boundaryPoint,'CORA:dimensionMismatch',I,[1;1],[0;1;0]);

% n-d arrays
lb = [];
lb(:,:,1,1) = [1 2; 3 5];
lb(:,:,1,2) = [0 -1; -2 3];
lb(:,:,1,3) = [1 1; -1 0];
lb(:,:,2,1) = [-3 2; 0 1];
ub = [];
ub(:,:,1,1) = [1.5 4; 4 10];
ub(:,:,1,2) = [1 2; 0 4];
ub(:,:,1,3) = [2 3; -0.5 2];
ub(:,:,2,1) = [-1 3; 0 2];
I = interval(lb,ub);
p = boundaryPoint(I,I.inf);
p_true = reshape([ 1.416667 4.000000 3.333333 8.333333 -2.500000 0.000000 2.833333 1.666667 0.500000 -1.333333 0.333333 4.000000 0.000000 0.000000 0.000000 0.000000 1.666667 -0.916667 2.166667 1.000000 0.000000 0.000000 0.000000 0.000000 ], [2,2,2,3]);
assert(all(withinTol(p, p_true,1e-6),"all"))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
