function res = test_levelSet_contains
% test_levelSet_contains - unit test function of levelSet containment
%
% Syntax:
%    res = test_levelSet_contains
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
% See also: -----

% Authors:       Maximilian Perschl
% Written:       08-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
% Define problem
% inclusion with linear levelSet
syms x y
eq = x + y - 1;
linear_ls = levelSet(eq,[x;y],'<=');
non_included_interval = interval([0.5;0],[1.5;1]);

% expect output false
assert(~contains(linear_ls,non_included_interval))

% different interval which should be included
included_interval = interval([0;0],[0.5;0.5]);

% expect output true
assert(contains(linear_ls,included_interval))

% inclusion with non-linear levelSet
eq = x^2 + y - 1;
nonlinear_ls = levelSet(eq,[x;y],'<=');

% expect output false
assert(~contains(nonlinear_ls,non_included_interval))

% expect output true
assert(contains(nonlinear_ls,included_interval))

% check for point inclusion close to border
% included
assert(contains(linear_ls,[0;1]))

% not included
assert(~contains(linear_ls,[0;1+1e-5],'exact',1e-6))

% multiple points
assert(all(contains(linear_ls,[0 0;1 1+1e-5],'exact',1e-6) == [1 0]));

% same for nonlinear set
% included
assert(contains(nonlinear_ls,[sqrt(2);-1],'exact',1e-5))

% not included
assert(~contains(nonlinear_ls,[sqrt(2);-1+1e-5],'exact',1e-6))

% multiple points
assert(all(contains(nonlinear_ls,[sqrt(2) sqrt(2);-1 -1+1e-5], ...
                                                'exact',1e-6) == [1 0]));

% inclusion with non-linear level set with multiple constraints
multiple_ls = levelSet([x - 1; x.^2 + y.^2 - 4],[x;y],'<=');
assert(contains(multiple_ls,[0;0]))
end

% ------------------------------ END OF CODE ------------------------------
