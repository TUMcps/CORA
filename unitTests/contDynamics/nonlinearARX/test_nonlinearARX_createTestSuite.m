function res = test_nonlinearARX_createTestSuite
% test_nonlinearARX_createTestSuite test for creating
%   trajectories for nonlinearARX
%
% Syntax:
%    res = test_nonlinearARX_createTestSuite
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%

% Authors:       Laura Luetzow
% Written:       28-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_k = 10;
n_m = 5;
n_s = 2;

% define system
dynamics = "NARX";
[sys, params.R0, params.U] = loadDynamics(dynamics);

% create trajectories and check dimensions
traj = createTestSuite(sys, params, n_k, n_m, n_s);
assert(length(traj) == n_m && traj(1).n_k == n_k && traj(1).n_s == n_s)

for inputCurve = ["rand", "randn", "bezier", "sigmoid","sinWave"]
    % specify options
    options.p_extr = 0.3;
    options.inputCurve = inputCurve;
    options.stateSet = 10*interval(-ones(sys.nrOfDims,1),...
        ones(sys.nrOfDims,1));

    % create trajectories and check dimensions
    traj = createTestSuite(sys, params, n_k, n_m, n_s, options);
    assert(length(traj) == n_m && traj(1).n_k == n_k && traj(1).n_s == n_s)

    % specify more options
    options2.inputCurve = inputCurve;
    options2.inputSet = interval([-100; 1; 0.5; 0; 1; 1],...
        [-1; 100; 1; 0.5; 2; 2]);

    % create trajectories and check dimensions
    traj = createTestSuite(sys, params, n_k, n_m, n_s, options);
    assert(length(traj) == n_m && traj(1).n_k == n_k && traj(1).n_s == n_s)
end

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
