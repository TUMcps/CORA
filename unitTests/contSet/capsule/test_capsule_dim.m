function res = test_capsule_dim
% test_capsule_dim - unit test function of dim
%
% Syntax:  
%    res = test_capsule_dim
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
% Written:      27-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Analytical -----------------------------------------------------------
% instantiate capsule
C = capsule([1;1],[1;1],0.5);

% compute enlarged capsule
C_dim = dim(C);

% true solution
true_dim = 2;

% compare results
res_analytical = C_dim == true_dim;
% -------------------------------------------------------------------------


% 2. Random ---------------------------------------------------------------
C = cell(10,1);
dimension = zeros(10,1);
for i=1:10
    dimension(i) = ceil(10*rand(1));
    center = rand(dimension(i),1);
    generator = rand(dimension(i),1);
    r = rand(1);
    C{i} = capsule(center,generator,r);
end

res_random = true;
for i=1:length(C)
    res_random = res_random && dim(C{i}) == dimension(i);
end
% -------------------------------------------------------------------------

% combine tests
res = res_analytical && res_random;

if res
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------