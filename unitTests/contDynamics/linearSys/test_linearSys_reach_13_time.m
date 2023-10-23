function res = test_linearSys_reach_13_time
% test_linearSys_reach_13_time - unit test function of linear reachability
%    analysis with shifted start time
%
% Syntax:
%    res = test_linearSys_reach_13_time
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       05-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

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

% model parameters
params.R0 = zonotope(10*ones(n,1),0.5*diag(ones(n,1)));
params.U = zonotope(zeros(3,1),0.05*diag(ones(3,1)));
params.tStart = 0.5;
params.tFinal = 1;

% call different algorithms
options.timeStep = 0.1;
options.taylorTerms = 4;
options.zonotopeOrder = 20;
algs = {'standard','wrapping-free','fromStart'};

% number of steps
steps = round((params.tFinal-params.tStart)/options.timeStep);

for i=1:length(algs)
    % set algorithm
    options.linAlg = algs{i};
    % reachability analysis
    R = reach(sys,params,options);

    % check if times are correct
    if ~withinTol(R.timePoint.time{1},params.tStart) ...
            || length(R.timePoint.time) ~= steps + 1 ...
            || ~withinTol(R.timePoint.time{end},params.tFinal)
        res = false;
    end

end

% adaptive algorithm
options_.linAlg = 'adaptive';
options_.error = 0.1;
% reachability analysis
R = reach(sys,params,options_);

% check if times are correct
if ~withinTol(R.timePoint.time{1},params.tStart) ...
        || ~withinTol(R.timePoint.time{end},params.tFinal)
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
