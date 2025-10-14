function completed = example_nonlinearSysDT_conform_02_gray
% example_nonlinearSysDT_conform_02_gray - example for conformance 
%   identification of nonlinear discrete-time systems to analyze the 
%   accuracy using simulated data
%
% Syntax:
%    example_nonlinearSysDT_conform_02_gray
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] L. Luetzow and M. Althoff, "Reachset-Conformant System
%        Identification," arXiv, 2025. 

% Authors:       Laura Luetzow
% Written:       14-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% User Specifications ----------------------------------------------------

% sytem dynamics: "lorenz","bicycle","bicycleHO","cstrDiscr","pedestrian"
dynamics = "lorenz"; 

% norm for evaluating the size of the reachable set: "interval","frob"
cost_norm = "interval"; 

% constraint type: "half", "gen"
constraints = "half";

% identification approach
methodsGray = ["graySeq","grayLS","graySim"];

% Set random number stream
rng(2)

% conformance identification
n_m = 3; % number of different input trajectories
n_s = 10; % number of sample trajectories per input trajectory
n_k = 3; % length of the identification trajectories

% validation
n_m_val = 2;
n_s_val = n_s;
n_k_val = n_k;

% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 1;
options_reach.tensorOrderOutput = 2;
options_reach.verbose = false;

% testSuite settings
options_testS.p_extr = 0.3;

% settings for plotting
plot_settings.plot_Yp = false; 
plot_settings.dims = [1 2]; 

% Load System Dynamics 
[sys, params_true.R0, params_true.U, p_true] = loadDynamics(dynamics,"rand");
params_true.tFinal = sys.dt * n_k - sys.dt;

% Create identification data 
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, ...
    n_s, options_testS);

%% Conformance Identification ---------------------------------------------

% Identification Options
options = options_reach;
options.cs.robustness = 1e-9;
options.cs.verbose = false;
options.cs.cost = cost_norm;
options.cs.constraints = constraints;
options.cs.p0 = 0.01*randn(size(p_true,1)+sys.nrOfInputs,1); % estimate parameter and center of U
options.cs.set_p = @(p,params) aux_set_p(p, params, dynamics);
options.cs.updateDeriv = false; % no recomputation of derivatives to decrease computation time

% Create struct for saving the identification results for each system
configs = cell(length(methodsGray) + 1,1);
configs{1}.sys = sys;
configs{1}.params = rmfield(params_true,'testSuite');
configs{1}.options = options_reach;
configs{1}.name = "true";


% Initial Estimates of the Disturbance Sets
c_R0 = center(params_true.R0);
c_U = center(params_true.U);

params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfDims) ones(sys.nrOfDims,1)]);
params_id_init.U = zonotope([c_U eye(sys.nrOfInputs) ...
    ones(sys.nrOfInputs,1)]);

% Identification ------------------------------------------------------

for i = 1:length(methodsGray)
    type = methodsGray(i);
    fprintf("Identification with method %s \n", type);
    timerVal = tic;
    [configs{i+1}.params, results] = conform(sys,params_id_init,options,type);
    Ts = toc(timerVal);
    configs{i+1}.sys = results.sys;
    configs{i+1}.options = options_reach;
    configs{i+1}.name = type;

    fprintf("Identification time: %.4f\n", Ts);
end

%% Validation and Visualization -------------------------------------------

% Sanity check: Compute Reachable Sets and Check Containment of the 
% Identification Test Cases
num_out = 0;
check_contain = 1;
methods = ["true" methodsGray];
for m=1:length(params_true.testSuite)
    [~, eval] = validateReach(params_true.testSuite(m), configs, check_contain);
    num_out = num_out + eval.num_out;
end
num_all = length(params_true.testSuite)*n_k_val*size(params_true.testSuite(1).y,3);
fprintf("IDENTIFICATION DATA: \n");
for i = 1:length(configs)
    p_contained = 100-num_out(i)/num_all*100;
    fprintf("%s: %.2f%% of the samples are contained in the reachable " + ...
        "set (must be 100%%!). \n", ...
        methods(i), p_contained);
end

% Create Validation Data
params_true.tFinal = sys.dt * n_k_val - sys.dt;
testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, ...
    n_s_val, options_testS);

% Compute Reachable Sets and Check Containment of the Validation Test Cases
num_out = 0;
check_contain = 1;
for m=1:length(testSuite_val)
    [~, eval] = validateReach(testSuite_val(m), configs, check_contain, plot_settings);
    num_out = num_out + eval.num_out;
end
num_all = length(testSuite_val)*n_k_val*size(testSuite_val(1).y,3);
fprintf("VALIDATION DATA: \n");
for i = 1:length(configs)
    p_contained = 100-num_out(i)/num_all*100;
    fprintf("%s: %.2f%% of the samples are contained in the " + ...
        "reachable set. \n", methods(i), p_contained);
end

% example completed
completed = true;

end


% Auxiliary functions -----------------------------------------------------

function [sys,params] = aux_set_p(p, params, dyn)
% return dynamics and uncertainty sets parameterized by p
[sys, ~, ~, p_t] = loadDynamics(dyn, 'diag', p);
c = p(length(p_t)+1:end);
params.U = zonotope(c, params.U.generators);
end

% ------------------------------ END OF CODE ------------------------------
