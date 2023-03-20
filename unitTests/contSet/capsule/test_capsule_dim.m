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
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      27-September-2019
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

%------------- END OF CODE --------------