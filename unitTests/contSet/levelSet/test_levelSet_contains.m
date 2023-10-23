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
if contains(linear_ls,non_included_interval)
    res = false;
end

% different interval which should be included
included_interval = interval([0;0],[0.5;0.5]);

% expect output true
if ~contains(linear_ls,included_interval)
    res = false;
end

% inclusion with non-linear levelSet
eq = x^2 + y - 1;
nonlinear_ls = levelSet(eq,[x;y],'<=');

% expect output false
if contains(nonlinear_ls,non_included_interval)
    res = false;
end

% expect output true
if ~contains(nonlinear_ls,included_interval)
    res = false;
end

% check for point inclusion close to border
% included
if ~contains(linear_ls,[0;1])
    res = false;
end
% not included
if contains(linear_ls,[0;1+1e-5],'exact',1e-6)
    res = false;
end
% same for nonlinear set
% included
if ~contains(nonlinear_ls,[sqrt(2);-1],'exact',1e-5)
    res = false;
end
% not included
if contains(nonlinear_ls,[sqrt(2);-1+1e-5],'exact',1e-6)
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
