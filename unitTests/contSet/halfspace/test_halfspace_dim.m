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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Analytical -----------------------------------------------------------
% instantiate capsule
h = halfspace([1;1],0.5);

% compute enlarged capsule
h_dim = dim(h);

% true solution
true_dim = 2;

% compare results
res_analytical = h_dim == true_dim;
% -------------------------------------------------------------------------


% 2. Random ---------------------------------------------------------------
h = cell(10,1);
dimension = zeros(10,1);
for i=1:10
    dimension(i) = ceil(10*rand(1));
    vector = rand(dimension(i),1);
    distance = 4*rand(1);
    h{i} = halfspace(vector,distance);
end

res_random = true;
for i=1:length(h)
    res_random = res_random && dim(h{i}) == dimension(i);
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