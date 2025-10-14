function res = test_trajectory_sequentialTrajectories()
% test_trajectory_sequentialTrajectories() - unit test function for creating
% sequential trajectories of a given length, which start 
% at every instant of the original trajectory (see Def. 7 in [1])
%
% Syntax:
%    test_trajectory_sequentialTrajectories()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 2022

% Authors:       Matthias Althoff, Laura Luetzow
% Written:       27-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% create traejctory
n_k = 30;
n_s = 2;
t = 1:n_k;
u = randn(2, n_k);
x = randn(3, n_k, n_s);
x(:,1,1) = x(:,1,2);
y = randn(2, n_k, n_s);
loc = randn(2, n_k);
a = randn(3, n_k, n_s);
n_k_new = 20;

traj = trajectory(u,x,y,t);
% try to create sequential trajectories with length n_k+5 (initial traj will be returned)
seqTraj = sequentialTrajectories(traj, n_k + 5);
assert(length(seqTraj) == 1);
assert(seqTraj.n_k == n_k);

% create sequential trajectories with length n_k_new
seqTraj = sequentialTrajectories(traj, n_k_new);
assert(length(seqTraj) == (n_k-n_k_new+1)*n_s);

% include loc and a
traj = trajectory(u,x,y,t,[],loc,a);
% create sequential trajectories
seqTraj = sequentialTrajectories(traj, n_k_new);
assert(length(seqTraj) == (n_k-n_k_new+1)*n_s);

% check correctness of each trajectory
for s = 0:n_s-1
    for i = 1:(n_k-n_k_new)
        m = (n_k-n_k_new+1)*s + i; % index of trajectory

        %% check if values starting from time step 2 equal first values of
        % next trajectory
        % check dimensions
        assert(seqTraj(m).n_k == n_k_new)
        assert(seqTraj(m).n_s == 1)

        % u
        assert(all(all(seqTraj(m).u(:,2:end,:) == seqTraj(m+1).u(:,1:end-1))));
        % x
        assert(all(all(seqTraj(m).x(:,2:end,:) == seqTraj(m+1).x(:,1:end-1,:))));
        % y
        assert(all(all(seqTraj(m).y(:,2:end,:) == seqTraj(m+1).y(:,1:end-1,:))));
        % u
        assert(all(all(seqTraj(m).u(:,2:end) == seqTraj(m+1).u(:,1:end-1))));
        % t
        assert(all(all(seqTraj(m).t(:,2:end) == seqTraj(m+1).t(:,1:end-1))));
        % loc
        assert(all(all(seqTraj(m).loc(:,2:end) == seqTraj(m+1).loc(:,1:end-1))));
        % a
        assert(all(all(seqTraj(m).a(:,2:end,:) == seqTraj(m+1).a(:,1:end-1,:))));
    end
end

% final result
res = true;

% ------------------------------ END OF CODE ------------------------------
