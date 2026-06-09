function res = test_neurNetContrSys_reach_04_allTypes
% test_neurNetContrSys_reach_04_allTypes - smoke tests for reach on
%    neurNetContrSys with linearSys, linearSysDT, linParamSys, and
%    nonlinearSysDT as the open-loop system
%
% Syntax:
%    res = test_neurNetContrSys_reach_04_allTypes
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Tobias Ladner
% Written:       20-February-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% shared neural network: 2 inputs -> 1 output
W1 = rand(10,2); b1 = rand(10,1);
W2 = rand(1,10); b2 = rand(1,1);
nn = neuralNetwork({ ...
    nnLinearLayer(W1, b1); nnReLULayer(); ...
    nnLinearLayer(W2, b2); nnReLULayer(); ...
});

% shared NN options
nn_options.bound_approx = true;
nn_options.poly_method = "singh";

% shared params (2-state systems, 1 NN-controlled input, short horizon)
dt = 0.05;
params.R0 = zonotope(interval([1; 0], [1.1; 0.1]));
params.U = zonotope(0);   % dummy external input (all inputs controlled by NN)
params.tFinal = 0.2;


% linearSys ----------------------------------------------------------

A = [-1 1; 0 -2]; B = [0; 1];
sysOL = linearSys(A, B);
sys = neurNetContrSys(sysOL, nn, dt);

options = struct();
options.maxError = Inf;
options.timeStep = dt;
options.taylorTerms = 4;
options.zonotopeOrder = 10;
options.nn = nn_options;

reach(sys, params, options);

% linearSysDT --------------------------------------------------------

A = [0.9 0.1; 0 0.8]; B = [0; 1];
sysOL = linearSysDT(A, B, dt);
sys = neurNetContrSys(sysOL, nn, dt);

options = struct();
options.maxError = Inf;
options.zonotopeOrder = 10;
options.nn = nn_options;

reach(sys, params, options);

% linParamSys --------------------------------------------------------

Ac = [-1 1; 0 -2]; Aw = [0.1 0; 0 0.1];
A = intervalMatrix(Ac, Aw); B = [0; 1];
sysOL = linParamSys(A, B, 'varParam');
sys = neurNetContrSys(sysOL, nn, dt);

options = struct();
options.maxError = Inf;
options.timeStep = dt;
options.taylorTerms = 4;
options.intermediateTerms = 4;
options.zonotopeOrder = 10;
options.nn = nn_options;

reach(sys, params, options);

% nonlinearSysDT -----------------------------------------------------

f = @(x,u) [x(1) + 0.05*x(2); x(2) + 0.05*u(1)];
sysOL = nonlinearSysDT(f, dt);
sys = neurNetContrSys(sysOL, nn, dt);

options = struct();
options.maxError = Inf;
options.zonotopeOrder = 10;
options.tensorOrder = 2;
options.errorOrder = 5;
options.lagrangeRem.zooMethods = 'interval';
options.nn = nn_options;

reach(sys, params, options);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
