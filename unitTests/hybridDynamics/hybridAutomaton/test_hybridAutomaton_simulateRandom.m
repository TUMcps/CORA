function res = test_hybridAutomaton_simulateRandom
% test_hybridAutomaton_simulateRandom - test function for simulateRandom
%
% Syntax:
%    res = test_hybridAutomaton_simulateRandom
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
% Written:       16-May-2023
% Last update:   15-October-2024 (MW, move to test_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% continuous dynamics 
A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];
linsys = linearSys('linearSys',A,B,c);

% rebound factor
alpha = -0.75;

% invariant set 
inv = polytope([-1,0],0);

% guard sets
guard = polytope([0,1],0,[1,0],0);

% reset function
reset = linearReset([0, 0; 0, alpha],zeros(2,1),zeros(2,1));

% self-transition
trans = transition(guard,reset,1);

% location object
loc = location('always',inv,trans,linsys); 

% hybrid automata
HA = hybridAutomaton(loc);

% model parameters
params.R0 = zonotope([[0.3;0],diag([0.025,0.025])]);
params.startLoc = 1;
params.finalLoc = 0;
params.tFinal = 0.5;

% settings for continuous reachability 
options.timeStep = 0.05;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% settings for hybrid systems
options.guardIntersect = 'polytope';
options.enclose = {'pca'};

% compute reachable set
R = reach(HA,params,options);

% simulate random trajectories
simOpt.points = 5;
simRes = simulateRandom(HA,params,simOpt);

% five trajectories
assert(length(simRes) == simOpt.points);

% correct start and end time
assert(all(arrayfun(...
    @(traj) traj.t{1}(1) == 0 && traj.t{end}(end) == params.tFinal,...
    simRes,'UniformOutput',true)));

% correct start point for each trajectory
assert(all(arrayfun(...
    @(traj) contains_(params.R0,traj.x{1}(1,:)','exact',eps),...
    simRes,'UniformOutput',true)));

% correct start location for each trajectory
assert(all(arrayfun(@(x) x.loc(1) == params.startLoc,...
    simRes,'UniformOutput',true)));

% time is strictly increasing
assert(all(arrayfun(@(x) all(diff(vertcat(x.t{:})) > -eps), ...
    simRes,'UniformOutput',true)));

% check if simulations are contained in reachable set
assert(contains(R,simRes));

% simulate edge case with a single input segment (caused a bug previously)
simOpt.nrConstInp = 1;
simulateRandom(HA,params,simOpt);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
