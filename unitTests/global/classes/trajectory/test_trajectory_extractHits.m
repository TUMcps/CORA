function res = test_trajectory_extractHits
% test_trajectory_extractHits - unit test function for extractHits
%
% Syntax:
%    res = test_trajectory_extractHits()
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

% TO-DO: adapt traj to hybrid systems where jump to different location
% creates new trajectory (with different length)

% empty trajectory
traj = trajectory();
[tHit,xHit,xHit_,iHit] = extractHits(traj);
assert(isempty(iHit) && isempty(tHit) && isempty(xHit) && isempty(xHit_));

% trajectory with trajectory (no specific location given)
n_k = 3;
t = 0.1*rand(1,n_k) + (1:n_k);
x = randn(2, n_k, 1);
traj = trajectory([],x,[],t);
[tHit,xHit,xHit_,iHit] = extractHits(traj);
assert(isempty(iHit) && isempty(tHit) && isempty(xHit) && isempty(xHit_));

% trajectory with one trajectory and given location
n_k = 6;
t = 0.1*rand(1,n_k) + (1:n_k);
x = randn(2, n_k, 1);
loc = [1 1 1 2 2 3];
traj = trajectory([],x,[],t,[],loc);

% no specific location
[tHit,xHit,xHit_,iHit] = extractHits(traj);
assert(all(iHit == [3 5]), all(withinTol(tHit, [t(3) t(5)], 1e-6)) ...
    && all(withinTol(xHit, [x(:,3,:) x(:,5,:)], 1e-6), 'all') ...
    && all(withinTol(xHit_, [x(:,4,:) x(:,6,:)], 1e-6), 'all'));
% extract specific location before jump
[tHit,xHit,xHit_,iHit] = extractHits(traj,3);
assert(isempty(iHit) && isempty(tHit) && isempty(xHit) && isempty(xHit_));
% extract specific location after jump
[tHit,xHit,xHit_,iHit] = extractHits(traj,[],2);
assert(iHit == 3, withinTol(tHit, t(3), 1e-6) ...
    && all(withinTol(xHit, x(:,3,:), 1e-6), 'all') ...
    && all(withinTol(xHit_, x(:,4,:), 1e-6), 'all'));

% trajectory with location vector
n_k = 6;
t = 0.1*rand(1,n_k) + (1:n_k);
x = randn(2, n_k, 1);
loc = [1 1 1 1 1 1; 1 1 1 2 2 2];
traj = trajectory([],x,[],t,[],loc);

% extract hits
[tHit,xHit,xHit_,iHit] = extractHits(traj,[1;1]);
assert(iHit == 3, withinTol(tHit, t(3), 1e-6) ...
    && all(withinTol(xHit, x(:,3,:), 1e-6), 'all') ...
    && all(withinTol(xHit_, x(:,4,:), 1e-6), 'all'));

% wrong size of location vector
[tHit,xHit,xHit_,iHit] = extractHits(traj,[1;1;1]);
assert(isempty(iHit) && isempty(tHit) && isempty(xHit) && isempty(xHit_));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
