function res = test_linearSys_simulateRandom_time
% test_linearSys_simulateRandom_time - unit test for random simulation
%    specifically checking whether start/end time and steps are correct
%
% Syntax:
%    res = test_linearSys_simulateRandom_time
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
% Written:       05-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
n = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];

% output matrix
C = [1 1 0 0;
     0 -0.5 0.5 0];

% initialize different linearSys-objects
sys = linearSys(A,B,zeros(4,1),C);

% model parameters and simulation settings
params.R0 = zonotope(10*ones(n,1),0.5*diag(ones(n,1)));
options.points = 5;

% only time horizon
params.tFinal = 1;
simRes = simulateRandom(sys,params,options);

% check if start and end time correct
res(end+1,1) = length(simRes) == options.points;
res(end+1,1) = all(arrayfun(@(z) z.t{1}(1) == 0,simRes,'UniformOutput',true));
res(end+1,1) = all(arrayfun(@(z) withinTol(z.t{1}(end),params.tFinal),simRes,'UniformOutput',true));
% check location
res(end+1,1) = all(arrayfun(@(z) isscalar(z.loc) && z.loc(1) == 0,simRes,'UniformOutput',true));

% set start time as well
params.tStart = 0.5;
params.tFinal = 1;
simRes = simulateRandom(sys,params,options);

% check if start and end time correct
res(end+1,1) = length(simRes) == options.points;
res(end+1,1) = all(arrayfun(@(z) withinTol(z.t{1}(1),params.tStart),simRes,'UniformOutput',true));
res(end+1,1) = all(arrayfun(@(z) withinTol(z.t{1}(end),params.tFinal),simRes,'UniformOutput',true));
% check location
res(end+1,1) = all(arrayfun(@(z) isscalar(z.loc) && z.loc(1) == 0,simRes,'UniformOutput',true));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
