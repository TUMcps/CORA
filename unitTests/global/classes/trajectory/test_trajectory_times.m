function res = test_trajectory_times
% test_trajectory_times - unit test function for times
%
% Syntax:
%    res = test_trajectory_times()
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

% Authors:       Tobias Ladner, Laura Luetzow
% Written:       08-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init trajectory object
n_k = 5;
n_s = 2;
t = 1:n_k;
x = randn(3, n_k, n_s);
x(:,1,1) = x(:,1,2);

traj = trajectory([],x,[],t);

% Column vector
A = [2; 3; -1];
traj_out = A .* traj;
assert(isequal(traj_out.x,A.*x));

% Scalar
A = 3;
traj_out = A .* traj;
assert(isequal(traj_out.x,A.*x));

A = 3;
traj_out = traj .* A;
assert(isequal(traj_out.x,A.*x));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
