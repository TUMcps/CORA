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

% Author:       Mark Wetzlinger
% Written:      28-August-2019
% Last update:  05-April-2023 (MW, unbounded intervals)
% Last revision:---

%------------- BEGIN CODE --------------

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

% check result (depends on ordering!)
res(end+1,1) = all(all(V == V_true));

% combine results
res = all(res);

%------------- END OF CODE --------------