function res = test_zonotope_boundaryPoint
% test_zonotope_boundaryPoint - unit test function of boundary point
%    computation
%
% Syntax:
%    res = test_zonotope_boundaryPoint
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
% Written:       16-April-2024
% Last update:   17-October-2024 (MW, integrate start point)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty zonotope
Z = zonotope.empty(2);
dir = [1;0];
x = boundaryPoint(Z,dir);
assert(all(size(x) == [2,0]));

% 1D zonotope
Z = zonotope(2,[1 -0.5, 1]);
% positive direction
x = boundaryPoint(Z,0.5);
x_true = 4.5;
assert(withinTol(x,x_true));
% negative direction
x = boundaryPoint(Z,-5);
x_true = -0.5;
assert(withinTol(x,x_true));
% different start point
x = boundaryPoint(Z,1,3);
x_true = 4.5;
assert(withinTol(x,x_true));

% 2D degenerate
Z = zonotope([1;-1],[2 4 -1; 1 2 -0.5]);
dir = [2;1];
x = boundaryPoint(Z,dir);
x_true = [8;2.5];
assert(all(withinTol(x,x_true)));
dir = [1;-2];
x = boundaryPoint(Z,dir);
x_true = Z.c;
assert(all(withinTol(x,x_true)));
% different start point
dir = [1;-2];
startPoint = [4.5;0.75];
x = boundaryPoint(Z,dir,startPoint);
assert(all(withinTol(x,startPoint)));

% 2D non-degenerate zonotope
Z = zonotope([1;-1],[-3 2 1; -1 0 3]);
dir = [1;1];
x = boundaryPoint(Z,dir);
x_true = [5;3];
assert(all(withinTol(x,x_true)));
% different start point
dir = [1;0];
startPoint = [-2;0];
x = boundaryPoint(Z,dir,startPoint);
x_true = [6;0];
assert(all(withinTol(x,x_true)));


% wrong calls
Z = zonotope([1;-1],[-3 2 1; -1 0 3]);
% ...all-zero 'direction'
assertThrowsAs(@boundaryPoint,'CORA:wrongValue',Z,[0;0]);
% ...start point not in the set
assertThrowsAs(@boundaryPoint,'CORA:wrongValue',Z,[1;1],[-500;100]);
% ...dimension mismatch
assertThrowsAs(@boundaryPoint,'CORA:dimensionMismatch',Z,[1;1;1],[-5;10]);
assertThrowsAs(@boundaryPoint,'CORA:dimensionMismatch',Z,[1;1],[0;1;0]);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
