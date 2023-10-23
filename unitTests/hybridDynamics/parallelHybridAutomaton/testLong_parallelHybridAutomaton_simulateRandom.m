function res = testLong_parallelHybridAutomaton_simulateRandom
% testLong_parallelHybridAutomaton_simulateRandom - test function for
%    random simulation
%
% Syntax:
%    res = testLong_parallelHybridAutomaton_simulateRandom
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

% system dynamics
PHA = roomHeatingParallel();

% simulation
params.startLoc = [1;1];
params.tFinal = 1;
params.R0 = zonotope([20.5;20.5],diag([0.1,0.1]));  
params.U = zonotope(4,0.01);

% settings for continuous reachability 
options.taylorTerms = 5; 
options.zonotopeOrder = 8; 
options.timeStep = 0.005; 
options.enclose = {'box','pca'}; 
options.guardIntersect = 'zonoGirard';

% reachability analysis
R = reach(PHA,params,options);

% simulation options
simOpt.points = 1;

% simulation
simRes = simulateRandom(PHA,params,simOpt);


% five resulting trajectories
res = length(simRes) == simOpt.points;
% reachable set has to contain trajectories
res(end+1,1) = contains(R,simRes);

res = all(res);

% ------------------------------ END OF CODE ------------------------------
