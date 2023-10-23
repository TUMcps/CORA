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
guard = conHyperplane([1 0],1);
reset = struct('A',zeros(2),'c',[-0.9;-0.9]);
trans(1) = transition(guard,reset,2);
guard = conHyperplane([0 1],1);
reset = struct('A',zeros(2),'c',[-0.5;0]);
trans(2) = transition(guard,reset,3);
flow = linearSys(zeros(2),0,[1;1]);
loc = location(inv,trans,flow);

% model parameters
params.loc = 1;
params.x0 = [-0.5;0];
params.tFinal = 2;

% simulate
[t,x,nextLoc,xJump] = simulate(loc,params);

% takes 1 second to hit guard
res = withinTol(t(end),1);
% hits guard at x = [0.5;1]
res(end+1,1) = compareMatrices(x(end,:)',[0.5;1]);
% all points must be in invariant
res(end+1,1) = all(contains(inv,x'));
% next location must be 3 (since second transition hit)
res(end+1,1) = nextLoc == 3;
% point after reset function
res(end+1,1) = compareMatrices(xJump,trans(2).reset.c);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
