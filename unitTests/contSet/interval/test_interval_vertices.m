function res = test_interval_vertices
% test_interval_vertices - unit test function of vertices
%
% Syntax:
%    res = test_interval_vertices
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       28-August-2019
% Last update:   05-April-2023 (MW, unbounded intervals)
%                28-April-2023 (MW, degenerate intervals)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-9;

% empty
I = interval.empty(2);
V = vertices(I);
res(end+1,1) = isnumeric(V) && isempty(V) && size(V,1) == 2;

% bounded
I = interval([-2; -4],[3; 1]);
V = vertices(I);
V_true = [-2 3  3 -2;
           1 1 -4 -4];
res(end+1,1) = compareMatrices(V,V_true,tol);

% unbounded
I = interval([-2; -4],[3; Inf]);
V = vertices(I);
V_true = [-2 -2   3 3;
          -4 Inf -4 Inf];
% check result (depends on ordering... cannot use compareMatrices here)
res(end+1,1) = all(all(V == V_true));

% degenerate
I = interval([-2; 0; 1],[5; 2; 1]);
V = vertices(I);
V_true = [-2 -2 5 5;
           0  2 0 2;
           1  1 1 1];
res(end+1,1) = compareMatrices(V,V_true);

% degenerate, point
lb = [1;4;-2;6];
I = interval(lb);
V = vertices(I);
res(end+1,1) = compareMatrices(V,lb);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
