function res = testLong_linearSys_reach_09_decomp
% testLong_linearSys_reach_09_decomp - unitTest for the decomposed
%    reachability algorithm from [1], compared to standard reach
%
% Syntax:
%    res = testLong_linearSys_reach_09_decomp
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References: 
%   [1] S. Bogomolov, M. Forets, G. Frehse, A. Podelski, C. Schilling
%       "Decomposing Reach Set Computations with Low-dimensional Sets and
%        High-Dimensional Matrices"

% Authors:       Mark Wetzlinger
% Written:       25-June-2019
% Last update:   14-August-2019
%                23-April-2020 (restructure params/options)
%                16-October-2024 (MW, changes due to proper decomp implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system dynamics
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;
fiveDimSys = linearSys('fiveDimSys',A,B);
C = [1 -1 0 0 0; -1 1 0 0 0];
fiveDimSys_outputs = linearSys('fiveDimSys',A,B,[],C);

% model parameters
params.tFinal = 1; 
params.R0 = zonotope(ones(fiveDimSys.dim,1),0.1*eye(fiveDimSys.dim));
params.U = 0.5*zonotope([1; 0; 0; 0.5; -0.5],diag([0.2, 0.5, 0.2, 0.5, 0.5]));

% reachability settings
options.timeStep = 0.1;
options.taylorTerms = 3;
options.zonotopeOrder = 50;
options.linAlg = 'decomp';
options.partition = [1,2; 3,4; 5,5];

% reachability: no C matrix -> output = states
R = reach(fiveDimSys, params, options);

% simulation
simOpt.points = 20;
traj = simulateRandom(fiveDimSys, params, simOpt);

% verify that simulation end points are contained in reachable set
finalSet = R.timePoint.set{end};
finalPoints = cell2mat(arrayfun(@(s) s.x(:,end),traj,'UniformOutput',false)');
assert(all(contains(finalSet,finalPoints,'exact',1e-8)));


% reachability with outputs
R = reach(fiveDimSys_outputs, params, options);

% simulation
simOpt.points = 20;
traj = simulateRandom(fiveDimSys_outputs, params, simOpt);

% verify that simulation end points are contained in reachable set
finalSet = R.timePoint.set{end};
finalPoints = cell2mat(arrayfun(@(s) s.y(:,end),traj,'UniformOutput',false)');
assert(all(contains(finalSet,finalPoints,'exact',1e-8)));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
