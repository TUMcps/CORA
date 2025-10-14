function res = test_trajectory_isemptyobject
% test_trajectory_isemptyobject - unit test function for isemptyobject
%
% Syntax:
%    res = test_trajectory_isemptyobject()
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

% Authors:       Mark Wetzlinger, Laura Luetzow
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty trajectory
traj = trajectory();
assert(isemptyobject(traj));

% trajectory with trajectory
t = [0 0.02 0.05];
x = [1 1 0.9; 1.1 0.8 1.2];
traj = trajectory([],x,[],t);
assert(~isemptyobject(traj));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
