function res = test_trajectory_find
% test_trajectory_find - unit test function for find
%
% Syntax:
%    res = test_trajectory_find()
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
% Written:       26-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% init trajectory with location
n_k = 5;
t = 1:n_k;
x = randn(3, n_k);
loc = [1 1 2 1 2];
traj = trajectory([],x,[],t,[],loc);

% find all points in location 1
traj_ = find(traj,'location',1);

% check correctness
assert(length(traj_) == 2 && all(traj_(1).t == [1 2]) ...
    && all(traj_(2).t == 4) && all(withinTol(traj_(1).x, x(:,1:2,:),1e-6),'all') ...
    && all(withinTol(traj_(2).x, x(:,4,:), 1e-6),'all') ...
    && all(traj_(1).loc == [1 1]) && all(traj_(2).loc == 1));

% same with location vector
loc = [1 1 2 1 2; 2 2 2 2 2];
traj = trajectory([],x,[],t,[],loc);

% find all points in location [1;2]
traj_ = find(traj,'location',[1;2]);

% check correctness
assert(length(traj_) == 2 && all(traj_(1).t == [1 2]) ...
    && all(traj_(2).t == 4) && all(withinTol(traj_(1).x, x(:,1:2,:),1e-6),'all') ...
    && all(withinTol(traj_(2).x, x(:,4,:), 1e-6),'all') ...
    && all(traj_(1).loc == [1 1; 2 2], 'all') && all(traj_(2).loc == [1; 2]));

% same with time
% find all points with time around 2
traj_ = find(traj,'time',2);

% check correctness
assert(length(traj_) == 1 && all(traj_.t == 2) ...
    && all(withinTol(traj_.x, x(:,2,:),1e-6),'all') ...
    && all(traj_.loc == [1; 2], 'all'));

% find all points with time in interval [1.5, 3.5]
traj_ = find(traj,'time',interval(1.5,3.5));

% check correctness
assert(length(traj_) == 1 && all(traj_.t == [2 3]) ...
    && all(withinTol(traj_.x, x(:,2:3,:),1e-6),'all') ...
    && all(traj_.loc == [1 2; 2 2], 'all'));

% no matching simulation should return isemptyobject
assert(isemptyobject(find(traj,'time',9)))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
