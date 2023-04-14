function res = test_halfspace_mtimes
% test_halfspace_mtimes - unit test function of mtimes
%
% Syntax:  
%    res = test_halfspace_mtimes
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
% Written:      16-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. empty case
% res_empty = true;
h = halfspace();
A = rand(2,2);
res_empty = isempty(A*h);

% 2. analytical test
% instantiate halfspace
c = [1;1]; d = 0.5;
h = halfspace(c,d);

% matrix
A = [2 3; 1 2];

% compute linear map
h_mapped = A*h;

% true result
A_inv = [2 -3; -1 2];
c_true = A_inv' * c;
h_true = halfspace(c_true,d);

% compare results
res_analytical = isequal(h_mapped,h_true);

% test singular matrix multiplication
try 
    A = [1 0; 0 0];
    h = A * h; % should throw an error
    res = false;
    return
end


% combine tests
res = res_empty && res_analytical;

%------------- END OF CODE --------------