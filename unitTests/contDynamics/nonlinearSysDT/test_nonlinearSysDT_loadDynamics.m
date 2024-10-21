function res = test_nonlinearSysDT_loadDynamics
% test_nonlinearSysDT_loadDynamics - unit test for loadDynamics function
%
% Syntax:
%    res = test_nonlinearSysDT_loadDynamics
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%

% Authors:       Laura Luetzow
% Written:       28-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng('default')

% Test different Uncertainty Set Types 
types = ["standard", "rand", "diag"];
for type = types
    [sys, params.R0, params.U, p_true] = loadDynamics("lorenz", type);

    % create test suite with no specification
    n_k = 10;
    n_m = 1;
    n_s = 1;
    T = createTestSuite(sys, params, n_k, n_m, n_s);
end

% Test different System Dynamics 
dynamics = ["cstrDiscr", "bicycle", "bicycleHO", "tank", "NARX", ...
    "Square", "pedestrian", "pedestrianARX"];
for dynamic = dynamics
    [sys, params.R0, params.U, p_true] = loadDynamics(dynamic);

    % create test suite with no specification
    n_k = 10;
    n_m = 1;
    n_s = 1;
    T = createTestSuite(sys, params, n_k, n_m, n_s);
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
