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
% Last update:  12-March-2021 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Empty case
C = capsule();

% compute dimension
C_dim = dim(C);

% true solution
true_dim = 0;

% compare results
res_empty = C_dim == true_dim;


% 2. Analytical
% instantiate capsule
C = capsule([1;1],[1;1],0.5);

% compute dimension
C_dim = dim(C);

% true solution
true_dim = 2;

% compare results
res_analytical = C_dim == true_dim;


% combine tests
res = res_empty && res_analytical;

if res
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------