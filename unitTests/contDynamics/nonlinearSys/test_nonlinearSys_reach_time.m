function res = test_nonlinearSys_reach_time
% test_nonlinearSys_reach_time - unit test function of nonlinear
%    reachability analysis with shifted start time
%
% Syntax:
%    res = test_nonlinearSys_reach_time
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

n = 3; m = 1;
f = @(x,u) [-x(2)*x(3); -x(1) + u(1); -x(2)*x(1)];
sys = nonlinearSys(f,n,m);

% model parameters
params.R0 = zonotope(2*ones(n,1),0.05*diag(ones(n,1)));
params.U = zonotope(0,0.01);
params.tStart = 0.15;
params.tFinal = 0.25;

% call different algorithms
options.timeStep = 0.01;
options.taylorTerms = 4;
options.zonotopeOrder = 20;
options.alg = 'lin';
options.tensorOrder = 2;

% number of steps
steps = round((params.tFinal-params.tStart)/options.timeStep);

% reachability analysis
R = reach(sys,params,options);

% check if times are correct
if ~withinTol(R.timePoint.time{1},params.tStart) ...
        || length(R.timePoint.time) ~= steps + 1 ...
        || ~withinTol(R.timePoint.time{end},params.tFinal)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
