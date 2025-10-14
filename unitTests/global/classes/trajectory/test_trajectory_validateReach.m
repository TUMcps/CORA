function res = test_trajectory_validateReach
% test_trajectory_validateReach- unit test function for evaluating and
%   plotting test cases
%
% Syntax:
%    test_trajectory_validateReach
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 2022

% Authors:       Laura Luetzow
% Written:       27-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% sytem dynamics: "lorenz","bicycle","bicycleHO","cstrDiscr"
for dynamics = ["pedestrian", "pedestrianARX", "lorenz", "NARX"]

    % Set random number stream
    rng(1)

    % conformance identification
    n_m = 3;
    n_s = 2;
    n_k = 6;

    % Reachability Settings
    options_reach.zonotopeOrder = 100;
    options_reach.tensorOrder = 2;
    options_reach.errorOrder = 5;
    options_reach.tensorOrderOutput = 2;

    % testSuite settings
    options_testS.p_extr = 0.8;

    % settings for plotting
    plot_settings.plot_Yp = false;
    plot_settings.dims = [1 2];

    % Load system dynamics and true uncertainty sets
    [sys, params_true.R0, params_true.U] = loadDynamics(dynamics,"rand");
    params_true.tFinal = sys.dt * n_k - sys.dt;

    % Create testSuite
    traj = createTestSuite(sys,params_true,n_k,n_m,n_s,options_testS);

    % Create struct for each system
    configs = cell(2,1);
    configs{1}.sys = sys;
    configs{1}.params = params_true;
    configs{1}.options = options_reach;
    configs{1}.name = "true";

    % Load system dynamics with different uncertainty sets
    [sys, params_true.R0, params_true.U] = loadDynamics(dynamics);
    configs{2}.sys = sys;
    configs{2}.params = params_true;
    configs{2}.options = options_reach;
    configs{2}.name = "different";

    % Initial Estimates of the Disturbance Sets
    c_R0 = center(params_true.R0);
    c_U = center(params_true.U);

    % Compute Reachable Sets
    [R, eval] = validateReach(traj(1), configs);

    % Compute Reachable Sets and Check Containment of the Test Cases
    num_out = 0;
    num_in = 0;
    check_contain = true;
    for m=1:length(traj)
        [~, eval] = validateReach(traj(m), configs, check_contain);
        num_out = num_out + eval.num_out;
        num_in = num_in + eval.num_in;
    end
    assert(num_out(1) == 0);
    assert(num_in(1) == n_m*n_k*n_s);
    if isa(sys, 'linearSysDT') || isa(sys, 'linearARX')
        [~, eval] = validateReach(traj, configs, check_contain);
    end

    % Compute Reachable Sets and Plot the Results
    check_contain = false;
    [R, eval] = validateReach(traj(1), configs, check_contain, plot_settings);
end

% example completed
res = true;
        
% ------------------------------ END OF CODE ------------------------------
