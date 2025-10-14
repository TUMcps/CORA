function res = test_trajectory_computeOutputDev
% test_trajectory_computeOutputDev - unit test function for computing the
%   output deviation y_a
%
% Syntax:
%    res = test_trajectory_computeOutputDev()
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

% Authors:       Laura Luetzow
% Written:       29-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_k = 10;
n_m = 2;
n_s = 2;

% define system
dynamics = "pedestrian";
[sys, params.R0, params.U] = loadDynamics(dynamics);

% create trajectories
traj = createTestSuite(sys, params, n_k, n_m, n_s);

% compute deviation
y_a = computeOutputDev(traj,sys);
assert(all(size(y_a, [1, 2]) == size(traj(1).y, [1, 2])));
assert(all(size(y_a, 3) == size(traj(1).y, 3)*n_m));

% define system
dynamics = "pedestrianARX";
[sys, params.R0, params.U] = loadDynamics(dynamics);

% create trajectories
traj = createTestSuite(sys, params, n_k, n_m, n_s);

% compute deviation
y_a = computeOutputDev(traj,sys);
assert(all(size(y_a, [1, 2]) == size(traj(1).y, [1, 2])));
assert(all(size(y_a, 3) == size(traj(1).y, 3)*n_m));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
