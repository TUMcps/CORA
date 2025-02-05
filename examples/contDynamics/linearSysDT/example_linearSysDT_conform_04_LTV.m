function completed = example_linearSysDT_conform_04_LTV
% example_linearSysDT_conform_04_LTV - example for conformance 
%   identification of nonlinear discrete-time systems to analyze the 
%   accuracy using simulated data
% 
%
% Syntax:
%    example_linearSysDT_conform_04_LTV
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% References:
%    [1]

% Authors:       Laura Luetzow
% Written:       31-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% User Specifications ----------------------------------------------------

dyn = "platoon"; % dynamics (choose from "platoon", "pedestrian")
n_k = 8; % number of time steps for identification
n_m = 20; % number of test cases for identification
n_s = 1; % number of samples per test case
n_n = 2; % number of vehicles in the platoon
n_m_val = 100; % number of test cases for validation
n_k_val = 20; % number of test cases for identification
constraints = {'gen','half',}; % type of containment constraints

% Set Random Number Stream
rng(2)

% Reachability Settings
options_reach.zonotopeOrder = 100;

% Conformance Settings
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.cost = 'interval';
options.cs.verbose = false;
options_testS.p_extr = 0.2; % probability of extreme points
options_testS.inputCurve = "rand"; % trajectory type for input

% Evaluation settings
check_contain = false;
plot_settings.dims = [1 2];
plot_settings.name = sprintf("Conformance Synthesis: %s", dyn);

% Parameters and System Dynamics
if dyn == "platoon"
    [sys, params_true.R0, params_true.U] = aux_load_platoon(n_n,...
        max(n_k,n_k_val));
else
    [sys, params_true.R0, params_true.U] = loadDynamics(dyn);
end
params_true.tFinal = sys.dt * n_k - sys.dt;

% Simulation
params_true.testSuite = createTestSuite(sys, ...
    params_true, n_k, n_m, n_s, options_testS);

% Initial Estimates of the Disturbance Sets
c_R0 = zeros(size(center(params_true.R0)));
c_U = zeros(size(center(params_true.U)));
params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfDims)]);
params_id_init.U = zonotope([c_U eye(sys.nrOfInputs)]);

%% Conformance Identification ---------------------------------------------
num_id = length(constraints);
name_id = cell(num_id,2);

% Struct for saving the identification results for each system
configs = cell(num_id+1,1);
configs{1}.sys = sys;
configs{1}.params = params_true;
configs{1}.options = options_reach;
configs{1}.name = "true";

for i_id = 1:num_id
    % run the identification
    options.cs.constraints = constraints{i_id};
    fprintf("Identification with %s-constraints, n_m=%d, " + ...
        "n_k=%d, n_x=%d\n",options.cs.constraints, n_m, n_k, 3*n_n);
    tic;
    [configs{i_id+1}.params, ~] = conform(sys,params_id_init,options);
    configs{i_id+1}.sys = sys;
    configs{i_id+1}.options = options_reach;
    configs{i_id+1}.name = options.cs.constraints;
    Ts=toc;
    fprintf("Identification time: %.4f\n", Ts);
end

%% Validation and Visualization -------------------------------------------

% Create Validation Data
if n_m_val ~= 0
    params_true.tFinal = sys.dt * n_k_val - sys.dt;
    testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, ...
        n_s, options);
end
% combine validation and trainings test cases
testSuite{1} = combineTestCases(params_true.testSuite{1}, testSuite_val{1});
plot_settings.s_val = size(params_true.testSuite{1}.y,3) + 1; 
    % (setting s_val leads to different color for validation test cases)

% run validation and plotting
validateReach(testSuite{1}, configs, check_contain, plot_settings);

% example completed
completed = true;

end


% Auxiliary functions -----------------------------------------------------

function [sys, R0, U] = aux_load_platoon(N_v,N_k)
% load the tank dynamics with the specified dimension

dt = 0.5;
N_u = N_v;
N_n = N_v*3;
sys = platoonN(dt,N_v,N_k);
c_R0 = randn(N_n,1); 
for i=0:N_v-1
    c_R0(i*N_v + 2) = 3*abs(c_R0(i*N_v + 2));
end
alpha_R0 = 2*rand(N_n,1);
c_U = randn(N_u,1);
alpha_U = rand(N_u,1);
R0 = zonotope([c_R0,diag(alpha_R0)]);
U = zonotope([c_U,diag(alpha_U)]);

end

% ------------------------------ END OF CODE ------------------------------
