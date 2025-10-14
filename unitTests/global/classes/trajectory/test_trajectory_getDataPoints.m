function res = test_trajectory_getDataPoints
% test_trajectory_getDataPoints - unit test function for transforming a
%   trajectory into data points
%
% Syntax:
%    res = test_trajectory_getDataPoints()
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
n_s = 1; % n_s > 1 is not implemented for getDataPoints
dim_x = 3;

t = (1:n_k) + 0.9*rand(1,n_k);
x = randn(dim_x, n_k, n_s);
traj = trajectory([],x,[],t);

% transform trajectory
points = getDataPoints(traj);
assert(size(points.x,1) == dim_x, size(points.x,2) == n_k-1)
assert(size(points.xNext,1) == dim_x, size(points.xNext,2) == n_k-1)
assert(all(withinTol(points.x(:,2:end), points.xNext(:,1:end-1), 1e-6), 'all'))
assert(size(points.dx,1) == dim_x, size(points.dx,2) == n_k-1)

% set compute_dx to true
points = getDataPoints(traj, true);

% set compute_dx to false
points = getDataPoints(traj, false);

% add input and output
dim_u = 2;
u = randn(dim_u, n_k);
y = randn(3, n_k, n_s);
traj = trajectory(u,x,y,t);

% convert trajectory
points = getDataPoints(traj);

% repeat for trajectory array
points = getDataPoints([traj; traj]);
assert(size(points.x,1) == dim_x, size(points.x,2) == 2*(n_k-1))
assert(size(points.xNext,1) == dim_x, size(points.xNext,2) == 2*(n_k-1))
assert(all(withinTol(points.x(:,2:n_k-1), points.xNext(:,1:n_k-2), 1e-6), 'all'))
assert(all(withinTol(points.x(:,n_k+1:end), points.xNext(:,n_k:end-1), 1e-6), 'all'))
assert(size(points.dx,1) == dim_x, size(points.dx,2) == 2*(n_k-1))
assert(size(points.u,1) == dim_u, size(points.u,2) == 2*(n_k-1))

% set compute_dx to true
points = getDataPoints([traj; traj], true);
assert(size(points.x,1) == dim_x, size(points.x,2) == 2*(n_k-1))
assert(size(points.xNext,1) == dim_x, size(points.xNext,2) == 2*(n_k-1))
assert(size(points.dx,1) == dim_x, size(points.dx,2) == 2*(n_k-1))
assert(size(points.u,1) == dim_u, size(points.u,2) == 2*(n_k-1))

% set compute_dx to false
points = getDataPoints([traj; traj], false);
assert(size(points.x,1) == dim_x, size(points.x,2) == 2*(n_k-1))
assert(size(points.xNext,1) == dim_x, size(points.xNext,2) == 2*(n_k-1))
assert(size(points.u,1) == dim_u, size(points.u,2) == 2*(n_k-1))

% result
res = true;

% ------------------------------ END OF CODE ------------------------------
