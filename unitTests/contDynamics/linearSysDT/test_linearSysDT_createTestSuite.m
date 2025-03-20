function res = test_linearSysDT_createTestSuite
% test_linearSysDT_createTestSuite - unit test for creating
%   test suites for given system dynamics
%
% Syntax:
%    res = test_linearSysDT_createTestSuite
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%

% Authors:       Laura Luetzow
% Written:       16-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System dimensions -------------------------------------------------------

nx = 3;
nu = 2;
ny = 2;

% Parameters --------------------------------------------------------------
params.R0 = zonotope(randn(nx,1),diag(randn(nx,1)));
params.U = zonotope(randn(nu,1),diag(randn(nu,1)));
rng('default')

% System Dynamics ---------------------------------------------------------

dt = 0.1;
A = randn(nx,nx);
B = randn(nx,nu);
C = randn(ny,nx);
D = randn(ny,nu);
sys = linearSysDT(A,B,[],C,D,dt);

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
options.stateSet = zonotope(3*ones(nx,1), 10*eye(nx));
options.inputSet = zonotope([2*sys.dt;0;n_k*sys.dt;8;8], diag([3*sys.dt, 1, 3*sys.dt, 10, 10]));
T = createTestSuite(sys, params, n_k, n_m, n_s, options);

% compute deviation between measurement trajectories and nominal output
% trajectories for test case m=1
T{1}=T{1}.compute_ya(sys);

res = true;
end

% ------------------------------ END OF CODE ------------------------------
