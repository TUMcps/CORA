function res = test_nonlinearSysDT_createTestSuite
% test_nonlinearSysDT_createTestSuite - unit test for creating
%   test suites for given system dynamics
%
% Syntax:
%    res = test_nonlinearSysDT_createTestSuite
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%

% Authors:       Laura Luetzow
% Written:       01-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------
params.R0 = zonotope([-0.15;-45],diag([0.005;3]));
params.U = zonotope(zeros(2,1),diag([0.1;2]));
rng('default')

% System Dynamics ---------------------------------------------------------

dt = 0.15;
fun = @(x,u) cstrDiscr(x,u,dt);
sys = nonlinearSysDT('stirredTankReactor',fun,0.015);

% create test suite with no specification
n_k = 100;
n_m = 5;
n_s = 20;
T = createTestSuite(sys, params, n_k, n_m, n_s);

% create test suite using the default parameter set
options.p_extr = 0.7;
options.inputFactor = [3 5];
for curve = ["rand", "randn", "bezier", "sinWave", "sigmoid"]
    options.inputCurve = curve;
    T = createTestSuite(sys, params, n_k, n_m, n_s, options);
end

% create test suite for given parameters
options.inputCurve = "bezier";
options.inputParameters = [[2*sys.dt;0;n_k*sys.dt;8;8], [1*sys.dt;-1;7*sys.dt;-8;-8]];
T = createTestSuite(sys, params, n_k, n_m, n_s, options);

% create test suite a given parameter set
options = rmfield(options, 'inputParameters');
options.inputCurve = "bezier";
options.stateSet = zonotope([3;3], 10*eye(2));
options.inputSet = zonotope([2*sys.dt;0;n_k*sys.dt;8;8], diag([3*sys.dt, 1, 3*sys.dt, 10, 10]));
T = createTestSuite(sys, params, n_k, n_m, n_s, options);


res = true;
end

% ------------------------------ END OF CODE ------------------------------
