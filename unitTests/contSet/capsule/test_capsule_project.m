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
cent = rand(dim,1);
gen = rand(dim,1);
rad = rand(1);
C = capsule(cent, gen, rad);

% project to lower dimension
projDim = [1 2 4];
C_proj = project(C, projDim);

% true projection
cent_true = cent(projDim);
gen_true  = gen(projDim);
% radius stays the same
C_true = capsule(cent_true, gen_true, rad);

% compare results
tol = 1e-9;
res_center = all(abs(center(C_proj) - center(C_true)) < tol);
res_generator = all(abs(C_proj.g - C_true.g) < tol);
res_radius = abs(radius(C_proj) - radius(C_true)) < tol;
res = res_center && res_generator && res_radius;

if res
    disp('test_project successful');
else
    disp('test_project failed');
end

%------------- END OF CODE --------------