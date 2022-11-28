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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      28-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate capsule
dim = 5;
c = rand(dim,1);
g = rand(dim,1);
r = rand(1);
C = capsule(c, g, r);

% project to lower dimension
projDim = [1 2 4];
C_proj1 = project(C, projDim);

% true projection
cent_true = c(projDim);
gen_true  = g(projDim);
% radius stays the same
C_true = capsule(cent_true, gen_true, r);

% projection using logical indices
projDim = [true true false true false];
C_proj2 = project(C, projDim);

% compare results
res = isequal(C_proj1,C_true) && isequal(C_proj2,C_true);

%------------- END OF CODE --------------