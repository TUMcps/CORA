function res = test_location_simulate
% test_location_simulate - test function for simulate
%
% Syntax:
%    res = test_location_simulate
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

% Authors:       Mark Wetzlinger
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init location with simple dynamics
inv = interval([-1;-1],[1;1]);
guard = polytope([],[],[1 0],1);
reset = linearReset(zeros(2),[],[-0.9;-0.9]);
trans(1) = transition(guard,reset,2);
guard = polytope([],[],[0 1],1);
reset = linearReset(zeros(2),[],[-0.5;0]);
trans(2) = transition(guard,reset,3);
flow = linearSys(zeros(2),[0;0],[1;1]);
loc = location(inv,trans,flow);

% model parameters
params.loc = 1;
params.x0 = [-0.5;0];
params.tFinal = 2;
% params.u = 0;

% simulate
[t,x,nextLoc,xJump] = simulate(loc,params);

% takes 1 second to hit guard
assert(withinTol(t(end),1));
% hits guard at x = [0.5;1]
assert(compareMatrices(x(:,end),[0.5;1]));
% all points must be in invariant
assert(all(contains(inv,x)));
% next location must be 3 (since second transition hit)
assert(nextLoc == 3);
% point after reset function
assert(compareMatrices(xJump,trans(2).reset.c));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
