function res = test_halfspace_dim
% test_halfspace_dim - unit test function of dim
%
% Syntax:  
%    res = test_halfspace_dim
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
% Last update:  16-March-2021 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

% 1. empty case
h = halfspace();
res_empty = true;
if dim(h) ~= 0
    res_empty = false;
end

% 2. analytical case
% instantiate halfspace
h = halfspace([1;1],0.5);

% compute dimension of halfspace
h_dim = dim(h);

% true solution
true_dim = 2;

% compare results
res_analytical = h_dim == true_dim;


% combine tests
res = res_empty && res_analytical;

if res
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------