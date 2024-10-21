function res = test_linearSys_reach_07_constInput
% test_linearSys_reach_07_constInput - unit test to check if constant
%    inputs c (cf. x' = Ax + Bu + c) are handled correctly; note: the
%    simulation results may be not with absolute certainty correct, but
%    should nonetheless remain inside the computed reachable sets
%
% Syntax:
%    test_linearSys_reach_07_constInput
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false whether reachable sets overapproximative
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       07-January-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% System dynamics ---------------------------------------------------------

A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];
linsys = linearSys('linearSys',A,B,c);

% Parameters --------------------------------------------------------------

params.R0 = zonotope([1;0],diag([0.05,0.05]));      % initial set
params.tFinal = 0.2;                                % final time
params.U = zonotope(0);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
options.taylorTerms = 10;
options.zonotopeOrder = 20;


% Reachability Analysis ---------------------------------------------------

R = reach(linsys,params,options);


% Simulation --------------------------------------------------------------

% number of initial points
simOpt.points = 5;
simRes = simulateRandom(linsys,params,simOpt); 


% Visualization -----------------------------------------------------------

% figure; hold on; box on;
% plot(R);
% plot(simRes);

% Numerical check ---------------------------------------------------------

% check containment of simulation results in reachable set
assert(contains(R,simRes));

% ------------------------------ END OF CODE ------------------------------
