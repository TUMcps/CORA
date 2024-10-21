function completed = example_nonlinearARX_conform_03_black
% example_nonlinearARX_conform_03_black - example for conformance 
%   identification of nonlinear discrete-time systems to analyze the 
%   accuracy using simulated data
%
% Syntax:
%    example_nonlinearARX_conform_03_black
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
dynamics = "Square"; 

% norm for evaluating the size of the reachable set: "interval","frob"
cost_norm = "interval"; 

% constraint type: "half", "gen"
constraints = "half";

% identification approach
methodsGray = ["blackGP","blackCGP"];

% Set random number stream
rng(2)

% conformance identification
n_m = 2; % number of different input trajectories
n_s = 50; % number of sample trajectories per input trajectory
n_k = 4; % length of the identification trajectories

% training and validation data
n_m_train = 100;
n_s_train = 10;
n_k_train = 4;
n_m_val = 5;
n_s_val = 10;
n_k_val = 4;

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
params_true.testSuite_train = createTestSuite(sys, params_true, ...
    n_k_train, n_m_train, n_s_train);
params_true.testSuite_val = createTestSuite(sys, params_true, ...
    n_k_val, n_m_val, n_s_val);

%% Conformance Identification ---------------------------------------------

% Identification Options
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.verbose = false;
options.cs.cost = cost_norm;
options.cs.constraints = constraints;

% Black-box approximation options
options.approx.gp_parallel = true;
options.approx.gp_pop_size = 50;
options.approx.gp_num_gen = 30;
options.approx.gp_func_names = {'times','plus', 'square'};
options.approx.gp_max_genes = 2;
options.approx.gp_max_depth = 2;
options.approx.gp_parallel = false;
options.approx.cgp_num_gen = 5;
options.approx.cgp_pop_size_base = 5;
options.approx.save_res = false;
options.approx.p = sys.n_p;

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
params_id_init.R0 = zonotope([c_R0]);
params_id_init.U = zonotope([c_U eye(size(c_U,1)) ones(size(c_U))]);

% Identification ------------------------------------------------------

for i = 1:length(methodsGray)
    type = methodsGray(i);
    fprintf("Identification with method %s \n", type);
    tic;
    [configs{i+1}.params, results] = conform(sys,params_id_init,options,type);
    Ts = toc;
    configs{i+1}.sys = results.sys;
    configs{i+1}.options = options_reach;
    configs{i+1}.name = type;

    fprintf("Identification time: %.4f\n", Ts);
end

%% Validation and Visualization -------------------------------------------

% Sanity check: Compute Reachable Sets and Check Containment of the 
% Identification Test Cases
num_out = 0;
num_in = 0;
check_contain = 1;
methods = ["true" methodsGray];
for m=1:length(params_true.testSuite)
    [~, eval] = validateReach(params_true.testSuite{m}, configs, check_contain);
    num_out = num_out + eval.num_out;
    num_in = num_in + eval.num_in;
end
num_all = length(params_true.testSuite)*n_k_val*size(params_true.testSuite{1}.y,3);
fprintf("IDENTIFICATION DATA: \n");
for i = 1:length(configs)
    p_contained = 100-(num_out(i)/(num_out(i)+num_in(i)))*100;
    fprintf("%s: %.2f%% of the samples are contained in the reachable " + ...
        "set (must be 100%%!). \n", ...
        methods(i), p_contained);
    fprintf("%s: %.2f%% of the samples were not valid (measurement " + ...
        "was nan or reachable set could not be computed). \n", ...
        methods(i), (num_all -(num_out(i)+num_in(i)))/num_all*100);
end

% Create Validation Data
params_true.tFinal = sys.dt * n_k_val - sys.dt;
testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, ...
    n_s_val, options_testS);

% Compute Reachable Sets and Check Containment of the Validation Test Cases
num_out = 0;
num_in = 0;
check_contain = 1;
for m=1:length(testSuite_val)
    [~, eval] = validateReach(testSuite_val{m}, configs, check_contain, plot_settings);
    num_out = num_out + eval.num_out;
    num_in = num_in + eval.num_in;
end
num_all = length(testSuite_val)*n_k_val*size(testSuite_val{1}.y,3);
fprintf("VALIDATION DATA: \n");
for i = 1:length(configs)
    p_contained = 100-(num_out(i)/(num_out(i)+num_in(i)))*100;
    fprintf("%s: %.2f%% of the samples are contained in the reachable " + ...
        "set. \n", ...
        methods(i), p_contained);
    fprintf("%s: %.2f%% of the samples were not valid. \n", ...
        methods(i), (num_all -(num_out(i)+num_in(i)))/num_all*100);
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
