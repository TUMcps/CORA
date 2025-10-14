function res = test_trajectory_uniformTimeStepSize
% test_trajectory_uniformTimeStepSize - unit test function for converting
%   trajectory to uniform time step size
%
% Syntax:
%    res = test_trajectory_uniformTimeStepSize()
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

% create trajectories
n_k = 5;
n_s = 2;

t = (1:n_k) + 0.9*rand(1,n_k);
x = randn(3, n_k, n_s);
x(:,1,1) = x(:,1,2);
traj = trajectory([],x,[],t);

% convert trajectory
dt = 1;
traj = uniformTimeStepSize(traj,dt);
assert(all(withinTol(0:n_k, traj.t,1e-6)))
assert(traj.dt == dt)
assert(traj.constant_dt)

% add input and output
u = randn(3, n_k);
y = randn(3, n_k, n_s);
traj = trajectory(u,x,y,t);

% convert trajectory
dt = 1;
traj = uniformTimeStepSize(traj,dt);
assert(all(withinTol(0:n_k, traj.t,1e-6)))
assert(traj.dt == dt)
assert(traj.constant_dt)

% repeat for trajectory array
traj = uniformTimeStepSize([traj; traj],dt);
assert(all(withinTol(0:n_k, traj(1).t,1e-6)) && all(withinTol(0:n_k, traj(2).t,1e-6)))
assert(traj(1).dt == dt && traj(2).dt == dt)
assert(traj(1).constant_dt && traj(2).constant_dt)

% don't add timestep size
traj = uniformTimeStepSize(traj);
assert(traj.constant_dt)

% repeat for trajectory array
traj = uniformTimeStepSize([traj; traj]);
assert(traj(1).constant_dt && traj(2).constant_dt)

% result
res = true;

% ------------------------------ END OF CODE ------------------------------
