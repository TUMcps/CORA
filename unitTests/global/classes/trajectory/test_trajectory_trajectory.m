function res = test_trajectory_trajectory
% test_trajectory_trajectory - unit test function for constructor
%
% Syntax:
%    res = test_trajectory_trajectory()
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

% empty trajectory
traj = trajectory();

for iter = 1:2
    if iter == 1
        n_s = 1;
    else
        n_s = 10;
    end

    % dimensions
    n_k = 5 + iter;
    dim_u = 2 + iter;
    dim_x = 4 + iter;
    dim_y = 3 + iter;
    dim_a = 3 + iter;

    % create random vectors
    u = randn(dim_u, n_k);
    x = randn(dim_x, n_k, n_s);
    x(:,1,:) = ones(dim_x,1,n_s); % same initial states
    y = randn(dim_y, n_k, n_s);
    a = randn(dim_a, n_k, n_s);
    t = 0.1*rand(1, n_k) + [1:n_k]; 
    dt = [];
    loc = ones(1,n_k);

    % correct instantiations using only the initial state
    traj = trajectory(u,x(:,1,1));
    traj = trajectory(u,x(:,1,1),y);
    traj = trajectory(u,x(:,1,1),y,t);
    traj = trajectory(u,x(:,1,1),y,t,dt,loc);
    traj = trajectory(u,x(:,1,1),y,t,dt,{},a);
    traj = trajectory(u,x(:,1,1),y,t,dt,{},[],'traj1');

    % correct instantiations with whole state vector
    traj = trajectory(u,x);
    traj = trajectory(u,x,y);
    traj = trajectory(u,x,y,t);
    traj = trajectory(u,x,y,t,dt,loc);
    traj = trajectory(u,x,y,t,dt,{},a);
    traj = trajectory(u,x,y,t,dt,{},[],'traj1');
end

% initialize sampling time
dt = 0.1;
t = [0:dt:(n_k-1)*dt];
traj = trajectory(u,x,y,t,dt);
traj = trajectory(u,x,y,[],dt);
assert(all(traj.t == t));
traj = trajectory(u,x,y,t);
assert(withinTol(traj.dt, dt, 1e-6));

% check wrong instantiations
% sampling time does not match the time vector
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x,y,t,2*dt);

% non-matching number of time steps
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u(:,3:end,:),x,y,t);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x(:,2:end,:),y,t);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x,y(:,2:end,:),t);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x,y,t(:,2:end,:));

% non-matching number of runs
u2 = randn(dim_u, n_k, n_s);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u2,x,y,t);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x(:,:,2:end),y,t);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x,y(:,:,2:end),t);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x,y,cat(3,t,t));

% non-matching initial states
x2 = randn(dim_x, n_k, n_s);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x2,y,t);
assertThrowsAs(@trajectory,'CORA:wrongInputInConstructor',u,x(:,1,:),y,t);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
