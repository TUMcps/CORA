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

% empty set
I_e = interval();
V_empty = vertices(I_e);
res = isnumeric(V_empty) && isempty(V_empty);

% create interval
lb = [-2; -4];
ub = [3; 1];
I = interval(lb, ub);

% compute vertices
V = vertices(I);

% true vertices
V_true = [-2 3  3 -2;
           1 1 -4 -4];

% check result
tol = 1e-9;
res(end+1,1) = compareMatrices(V,V_true,tol);

% unbounded interval
lb = [-2; -4];
ub = [3; Inf];
I = interval(lb, ub);

% compute vertices
V = vertices(I);

% true vertices
V_true = [-2 -2   3 3;
          -4 Inf -4 Inf];

% check result (depends on ordering... cannot use compareMatrices here)
res(end+1,1) = all(all(V == V_true));


% degenerate interval
lb = [-2; 0; 1];
ub = [5; 2; 1];
I = interval(lb,ub);

% compute vertices
V = vertices(I);

% true vertices
V_true = [-2 -2 5 5;
           0  2 0 2;
           1  1 1 1];

% check result
res(end+1,1) = compareMatrices(V,V_true);


% interval is just a point
lb = [1;4;-2;6];
I = interval(lb);

% compute vertices
V = vertices(I);

% check result
res(end+1,1) = compareMatrices(V,lb);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
