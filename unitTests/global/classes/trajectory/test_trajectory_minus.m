function res = test_trajectory_minus
% test_trajectory_minus - unit test function for minus
%
% Syntax:
%    res = test_trajectory_minus()
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
dim_x = 3;
n_k = 100;

x = rand(dim_x,n_k,2);
x(:,1,2) = x(:,1,1); % same initial state
t = 1:n_k;

traj = trajectory([],x,[],t);

% constant bias
shift = 2;
traj_out = traj - shift;
assert(isequal(traj_out.x,x-shift));

shift = 2;
traj_out = shift - traj;
assert(isequal(traj_out.x,shift-x));

% same shift for every trajectory and time step
shift = [0.1;0.2;0.3];
traj_out = traj - shift;
assert(isequal(traj_out.x,x-shift));

% same shift for every trajectory
shift = rand(dim_x,n_k);
traj_out = shift - traj;
assert(isequal(traj_out.x,shift-x));

% custom shift for every state
shift = rand(dim_x,n_k,2);
traj_out = shift - traj;
assert(isequal(traj_out.x,shift-x));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
