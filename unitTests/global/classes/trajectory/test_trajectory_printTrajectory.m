function res = test_trajectory_printTrajectory
% test_trajectory_printTrajectory - unit test function of printSpec
%
% Syntax:
%    res = test_trajectory_printTrajectory
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

% test print of a simple system
n_k = 5;
n_s = 2;
t = 1:n_k;
x = randn(3, n_k, n_s);
x(:,1,1) = x(:,1,2);
traj = trajectory([],x,[],t);

printTrajectory(traj);
printTrajectory(traj,'high');
printTrajectory(traj,'high',true);
printTrajectory(traj,'high',false);

% test fid
filename = 'test.txt';
printTrajectory(filename,traj,'high',true);
traj_copy = eval(fileread(filename));
% assert(isequal(traj,traj_copy)); % not implemented
delete(filename)
assert(all(withinTol(traj.x,traj_copy.x,1e-6),'all'))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
