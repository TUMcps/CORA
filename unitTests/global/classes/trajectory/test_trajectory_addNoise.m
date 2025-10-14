function res = test_trajectory_addNoise
% test_trajectory_addNoise - unit test function for adding noise to the
%       state of a trajectory
%
% Syntax:
%    res = test_trajectory_addNoise()
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
% Written:       25-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% dimensions
n_k = 5;
n_s = 5;
dim_u = 2;
dim_x = 4;
dim_y = 3;
dim_a = 3;

traj(5,1) = trajectory();
for m = 1:5

    % create random vectors
    u = randn(dim_u, n_k, 1);
    x = randn(dim_x, n_k, n_s);
    x(:,1,:) = ones(dim_x, 1, n_s);
    y = randn(dim_y, n_k, n_s);
    a = randn(dim_a, n_k, n_s);
    t = [1:n_k];
    dt = [];
    loc = ones(1, n_k);

    % instantiation
    traj(m) = trajectory(u,x,y,t,dt,loc,a,'traj1');
end
% add noise to the trajectories and chack if initial state is consistent
trajectoriesNoise = addNoise(traj,0.2, 'x');
assert(all(abs(diff(trajectoriesNoise(1).x(:,1,:),[],3)) < 1e-6,'all'))
trajectoriesNoise = addNoise(traj,1, 'xy','uniform');
assert(all(abs(diff(trajectoriesNoise(2).x(:,1,:),[],3)) < 1e-6,'all'))
trajectoriesNoise = addNoise(traj,0.5, 'xyu','gaussian');
assert(all(abs(diff(trajectoriesNoise(5).x(:,1,:),[],3)) < 1e-6,'all'))
trajectoriesNoise = addNoise(traj,0.5, 'yua','gaussian');
assert(all(abs(diff(trajectoriesNoise(5).x(:,1,:),[],3)) < 1e-6,'all'))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
