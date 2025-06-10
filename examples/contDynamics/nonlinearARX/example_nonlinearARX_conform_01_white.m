function completed = example_nonlinearARX_conform_01_white
% example_nonlinearARX_conform_01_white - example for conformance 
%   identification of nonlinear ARX systems to analyze the accuracy using
%   simulated data
%
% Syntax:
%    example_nonlinearARX_conform_01_white
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] ???

% Authors:       Laura Luetzow
% Written:       03-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% User Specifications ----------------------------------------------------

% sytem dynamics: "NARX","NARX2"
dynamics = "NARX"; 

% norm for evaluating the size of the reachable set: "interval","frob"
cost_norm = "interval"; 

% constraint type: "half", "gen"
constraints = "half";

% Set random number stream
rng(1)

% conformance identification
n_m = 30; % number of different input trajectories
n_s = 100; % number of sample trajectories per input trajectory
n_k = 4; % length of the identification trajectories

% validation
n_m_val = 2;
n_s_val = n_s;
n_k_val = n_k;

% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 1;
options_reach.tensorOrderOutput = 2;

% testSuite settings
options_testS.p_extr = 0.3;

% settings for plotting
plot_settings.plot_Yp = false; 
plot_settings.dims = [1 2]; 

% Load System Dynamics 
[sys, params_true.R0, params_true.U] = loadDynamics(dynamics,"rand");
params_true.tFinal = sys.dt * n_k - sys.dt;

% Create identification data 
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, ...
    n_s, options_testS);

%% Conformance Identification ---------------------------------------------

% Identification Options
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.verbose = false;
options.cs.cost = cost_norm;
options.cs.constraints = constraints;

% Create struct for saving the identification results for each system
configs = cell(2,1);
configs{1}.sys = sys;
configs{1}.params = rmfield(params_true,'testSuite');
configs{1}.options = options_reach;
configs{1}.name = "true";

% Initial Estimates of the Disturbance Sets
c_R0 = center(params_true.R0);
c_U = center(params_true.U);

params_id_init = params_true;
params_id_init.R0 = zonotope(c_R0);
params_id_init.U = zonotope([c_U eye(size(c_U,1)) ones(size(c_U))]);

% Identification ------------------------------------------------------
fprintf("Identification with %s-cost, %s-constraints \n", ...
    options.cs.cost, options.cs.constraints);
timerVal = tic;
[configs{2}.params] = conform(sys,params_id_init,options,"white");
Ts = toc(timerVal);
configs{2}.sys = sys;
configs{2}.options = options_reach;
configs{2}.name = options.cs.constraints;

fprintf("Identification time: %.4f\n", Ts);

%% Validation and Visualization -------------------------------------------

% Sanity check: Compute Reachable Sets and Check Containment of the 
% Identification Test Cases
num_out = 0;
check_contain = 1;
for m=1:length(params_true.testSuite)
    [~, eval] = validateReach(params_true.testSuite{m}, configs, check_contain);
    num_out = num_out + eval.num_out;
end
num_all = length(params_true.testSuite)*n_k_val*size(params_true.testSuite{1}.y,3);
p_contained = 100-num_out(2)/num_all*100;
fprintf("%.2f%% of the samples are contained in the reachable " + ...
    "set for the identification test cases (must be 100%%!). \n", ...
    p_contained);

% Create Validation Data
params_true.tFinal = sys.dt * n_k_val - sys.dt;
testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, ...
    n_s_val, options_testS);

% Compute Reachable Sets and Check Containment of the Validation Test Cases
num_out = 0;
check_contain = 1;
for m=1:length(testSuite_val)
    [~, eval] = validateReach(testSuite_val{m}, configs, check_contain, plot_settings);
    num_out = num_out + eval.num_out;
end
num_all = length(testSuite_val)*n_k_val*size(testSuite_val{1}.y,3);
p_contained = 100-num_out(2)/num_all*100;
fprintf("%.2f%% of the samples are contained in the " + ...
    "reachable set for the validation test cases. \n", p_contained);

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
