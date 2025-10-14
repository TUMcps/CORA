function res = test_trajectory_mtimes
% test_trajectory_mtimes - unit test function for mtimes
%
% Syntax:
%    res = test_trajectory_mtimes()
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
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init trajectory object
n_k = 5;
n_s = 2;
t = 1:n_k;
x = randn(3, n_k, n_s);
x(:,1,1) = x(:,1,2); % same initial state

traj = trajectory([],x,[],t);

% simple cases
A = [2 3 -1; 1 2 4];
traj_out = A * traj;
assert(all(withinTol(traj_out.x(:,:,1),A*x(:,:,1),1e-9),'all') ...
    && all(withinTol(traj_out.x(:,:,2),A*x(:,:,2),1e-9),'all'));

A = 3;
traj_out = A * traj;
assert(all(withinTol(traj_out.x(:,:,1),A*x(:,:,1),1e-9),'all') ...
    && all(withinTol(traj_out.x(:,:,2),A*x(:,:,2),1e-9),'all'));

A = 3;
traj_out = traj * A;
assert(all(withinTol(traj_out.x(:,:,1),A*x(:,:,1),1e-9),'all') ...
    && all(withinTol(traj_out.x(:,:,2),A*x(:,:,2),1e-9),'all'));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
