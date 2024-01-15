function res = test_capsule_project
% test_capsule_project - unit test function of project
%
% Syntax:
%    res = test_capsule_project
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
% Written:       28-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% instantiate capsule
c = [1; 0; -1; 2; 1];
g = [0; 0; -1; 1; 2];
r = 0.75;
C = capsule(c, g, r);

% project to lower dimension
projDim = [1 2 4];
C_proj1 = project(C, projDim);
% projection using logical indices
projDim = [true true false true false];
C_proj2 = project(C, projDim);

% true projection (same radius)
c_true = c(projDim);
g_true  = g(projDim);
C_true = capsule(c_true, g_true, r);

res(end+1,1) = isequal(C_proj1,C_true);
res(end+1,1) = isequal(C_proj2,C_true);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
