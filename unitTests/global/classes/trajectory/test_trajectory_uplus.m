function res = test_trajectory_uplus
% test_trajectory_uplus - unit test function for uplus
%
% Syntax:
%    res = test_trajectory_uplus()
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

% simple cases
traj_out = +traj;
assert(isequal(traj_out.x,x));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
