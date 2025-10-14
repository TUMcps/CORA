function res = test_trajectory_add
% test_trajectory_add - unit test function for adding two trajectory
%       objects together
%
% Syntax:
%    res = test_trajectory_add()
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
traj = add(traj,traj);

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

    u = randn(dim_u, n_k, 1);
    t = 0.1*rand(1, n_k, 1) + [1:n_k]; 

    % create random vectors for traj1
    x1 = randn(dim_x, n_k, n_s);
    x1(:, 1, :) = ones(dim_x, 1, n_s);
    y1 = randn(dim_y, n_k, n_s);
    a1 = randn(dim_a, n_k, n_s);

    % create random vectors for traj2
    x2 = randn(dim_x, n_k, n_s);
    x2(:, 1, :) = x1(:, 1, :);
    y2 = randn(dim_y, n_k, n_s);
    a2 = randn(dim_a, n_k, n_s);

    dt = [];
    loc = ones(1,n_k);

    % correct instantiations using only the initial state
    traj1 = trajectory(u,x1(:,1,1));
    traj2 = trajectory(u,x2(:,1,1));
    traj = add(traj1,traj2);

    traj1 = trajectory(u,x1(:,1,1),y1,t,dt,loc);
    traj2 = trajectory(u,x2(:,1,1),y2,t,dt,loc);
    traj = add(traj1,traj2);

    traj1 = trajectory(u,x1(:,1,1),y1,t,dt,[],a1,'traj1');
    traj2 = trajectory(u,x2(:,1,1),y2,t,dt,[],a2,'traj2');
    traj = add(traj1,traj2);

    % correct instantiations with whole state vector
    traj1 = trajectory(u,x1,y1,t);
    traj2 = trajectory(u,x2,y2,t);
    traj = add(traj1,traj2);
end

% check wrong instantiations
% non-matching number of time steps
traj1 = trajectory(u(:,1:end-1,:),x1(:,1:end-1,:),y1(:,1:end-1,:),t(:,1:end-1,:));
traj2 = trajectory(u,x2,y2,t);
assertThrowsAs(@add,'CORA:specialError',traj1,traj2);

% non-matching of non-empty vectors
traj1 = trajectory(u,[],y1,t);
traj2 = trajectory(u,x2,y2,t);
assertThrowsAs(@add,'CORA:specialError',traj1,traj2);

% not the same input vectors
u2 = randn(dim_u, n_k, 1);
traj1 = trajectory(u,x1,y1,t);
traj2 = trajectory(u2,x2,y2,t);
assertThrowsAs(@add,'CORA:specialError',traj1,traj2);

% not the same initial states
traj1 = trajectory(u,x1,y1,t);
traj2 = trajectory(u,0.1* x2,y2,t);
assertThrowsAs(@add,'CORA:specialError',traj1,traj2);

% not the same time vectors
t2 = 0.01*rand(1, n_k, 1) + [1:n_k]; 
traj1 = trajectory(u,x1,y1,t);
traj2 = trajectory(u,x2,y2,t2);
assertThrowsAs(@add,'CORA:specialError',traj1,traj2);

% not the same location vectors
traj1 = trajectory(u,x1,y1,t,dt,loc);
traj2 = trajectory(u,x2,y2,t,dt,2*loc);
assertThrowsAs(@add,'CORA:specialError',traj1,traj2);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
