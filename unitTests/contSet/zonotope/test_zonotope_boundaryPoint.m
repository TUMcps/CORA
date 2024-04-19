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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty zonotope
Z = zonotope.empty(2);
dir = [1;0];
x = boundaryPoint(Z,dir);
res(end+1,1) = all(size(x) == [2,0]);

% 1D zonotope
Z = zonotope(2,[1 -0.5, 1]);
dir = 0.5;
x = boundaryPoint(Z,dir);
x_true = 4.5;
res(end+1,1) = withinTol(x,x_true);
dir = -5;
x = boundaryPoint(Z,dir);
x_true = -0.5;
res(end+1,1) = withinTol(x,x_true);

% 2D degenerate
Z = zonotope([1;-1],[2 4 -1; 1 2 -0.5]);
dir = [2;1];
x = boundaryPoint(Z,dir);
x_true = [8;2.5];
res(end+1,1) = all(withinTol(x,x_true));
dir = [1;-2];
x = boundaryPoint(Z,dir);
x_true = Z.c;
res(end+1,1) = all(withinTol(x,x_true));

% 2D non-degenerate zonotope
Z = zonotope([1;-1],[-3 2 1; -1 0 3]);
dir = [1;1];
x = boundaryPoint(Z,dir);
x_true = [5;3];
res(end+1,1) = all(withinTol(x,x_true));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
