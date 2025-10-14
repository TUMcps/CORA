function res = testLong_linearSys_reach_11_inputTrajectory
% testLong_linearSys_reach_11_inputTrajectory - unit test for linear
%    reachability analysis, different algorithms compute the reachable
%    sets given time-varying input signal; afterwards, the trajectories
%    are checked for containment within the reachable sets
%
% Syntax:
%    res = testLong_linearSys_reach_11_inputTrajectory
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
% Written:       16-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------

A = [-0.7 -2; 2 -0.7];
B = 1;

sys = linearSys('sys',A,B);


% Parameters --------------------------------------------------------------

dim_x = length(A);

params.tFinal = 5;
params.R0 = zonotope([[10; 5],0.5*eye(dim_x)]);
params.U = zonotope([zeros(dim_x,1),0.25*eye(dim_x)]);
params.u = 50*[0.01 1; -0.02 -0.5];


% Simulation --------------------------------------------------------------

simOpt.points = 2;
simOpt.fracVert = 1;
simOpt.fracInpVert = 1;

traj = simulateRandom(sys, params, simOpt);


% Reachability Analysis ---------------------------------------------------

options.linAlg = 'adaptive';
options.error = 0.05;
Radaptive = reach(sys,params,options);

% reachability settings for non-adaptive algorithms
options = rmfield(options,'error');
options.timeStep = 0.05;
options.taylorTerms = 5;
options.zonotopeOrder = 100;

% correct dimension for u
params.u = repelem(params.u,1,...
    params.tFinal/options.timeStep/size(params.u,2));

% standard algorithm
options.linAlg = 'standard';
Rstandard = reach(sys,params,options);

% wrapping free
options.linAlg = 'wrapping-free';
Rwrappingfree = reach(sys,params,options);


% Verification ------------------------------------------------------------

assert(contains(Radaptive,traj));
assert(contains(Rstandard,traj));
assert(contains(Rwrappingfree,traj));

res = true;

% ------------------------------ END OF CODE ------------------------------
