function res = test_nonlinearSysDT_conform_01_constraints
% test_nonlinearSysDT_conform_01_constraints - unit test for comparing
%        reachset-conformant identification with different constraints
%
% Syntax:
%    res = test_nonlinearSysDT_conform_01_constraints
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       06-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% set random number stream
rng('default')

cost_norm = "interval"; 
constraints = {"half", "gen"};
n_m = 2;
n_s = 30;
n_k = 4;

% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 5;
options_reach.tensorOrderOutput = 2;

% Load System Dynamics 
[sys, params_true.R0, params_true.U] = loadDynamics('lorenz');
params_true.tFinal = sys.dt * n_k - sys.dt;

% Create identification data 
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s);

%% Conformance Identification ---------------------------------------------

% General Options
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.verbose = false;
options.cs.cost = cost_norm;

num_id = length(constraints);
params_id = cell(num_id,1);
results = cell(num_id,1);

for i_id = 1:num_id
    options.cs.constraints = constraints{i_id};

    % Initial Estimates of the Disturbance Sets
    c_R0 = zeros(size(center(params_true.R0)));
    c_U = zeros(size(center(params_true.U)));
    params_id_init = params_true;
    params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfDims) ones(sys.nrOfDims,1)]);
    params_id_init.U = zonotope([c_U eye(sys.nrOfInputs) ...
        ones(sys.nrOfInputs,1)]);

    % Identification ------------------------------------------------------
    [params_id{i_id}, results{i_id}] = conform(sys,params_id_init,options);   
end

G_X01 = params_id{1}.R0.G;
G_X02 = params_id{2}.R0.G;
G_U1 = params_id{1}.U.G;
G_U2 = params_id{2}.U.G;

tolX = max(1e-1*max([G_X01 G_X02], [],'all'),1e-6); 
tolU = max(1e-1*max([G_U1 G_U2], [],'all'),1e-6); 

% should lead to the same cost value and approximately the same uncertainty
% sets
assert(abs(results{1}.fval - results{2}.fval)<= 1e-5)
assert(all(abs(G_X01-G_X02) <= tolX, 'all'))
assert(all(abs(G_U1-G_U2) <= tolU, 'all'))

res = true;

end

% ------------------------------ END OF CODE ------------------------------
