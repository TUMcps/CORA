function res = testLong_nonlinearARX_conform_03_black
% testLong_nonlinearARX_conform_03_black - unit test for testing the
%        reachset-conformant black-box identification approaches
%
% Syntax:
%    res = testLong_nonlinearARX_conform_03_black
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       31-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% set random number stream
rng('default')

cost_norm = "interval"; 
constraints = "half";
id_methods = ["blackGP","blackCGP"];   
n_m = 2;
n_s = 50;
n_k = 4;
n_m_train = 100;
n_s_train = 10;
n_k_train = 4;
n_m_val = 5;
n_s_val = 10;
n_k_val = 4;


% Load System Dynamics 
[sys, params_true.R0, params_true.U] = loadDynamics('NARX');
params_true.tFinal = sys.dt * n_k - sys.dt;

% Create identification data 
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s);
params_true.testSuite_train = createTestSuite(sys, params_true, ...
    n_k_train, n_m_train, n_s_train);
params_true.testSuite_val = createTestSuite(sys, params_true, ...
    n_k_val, n_m_val, n_s_val);

% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 5;
options_reach.tensorOrderOutput = 2;

% Conformance Options
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.verbose = false;
options.cs.cost = cost_norm;
options.cs.constraints = constraints;

% Black-box approximation options
options.approx.nn_lr = 5e-3;
options.approx.gp_parallel = true;
options.approx.gp_pop_size = 50;
options.approx.gp_num_gen = 20;
options.approx.cgp_num_gen = 2;
options.approx.gp_parallel = false;
options.approx.cgp_pop_size_base = 2;
options.approx.save_res = false;
options.approx.p = sys.n_p;
% options.approx.nn_lrDropPeriod = 6;
% options.approx.nn_lrDropFactor = 0.1;
% options.approx.nn_act = "tanh";
% options.approx.nn_bs = 256;
% options.approx.nn_neurons = [16 16];
% options.approx.cgp_n_m_conf = 5;
% options.approx.folder_data = folder_data;


%% Conformance Identification ---------------------------------------------

% Initial Estimates of the Disturbance Sets
c_R0 = zeros(size(center(params_true.R0)));
c_U = zeros(size(center(params_true.U)));
params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0 eye(sys.n_p*sys.nrOfOutputs) ones(sys.n_p*sys.nrOfOutputs,1)]);
params_id_init.U = zonotope([c_U eye(sys.nrOfInputs) ...
    ones(sys.nrOfInputs,1)]);

for i_id = 1:length(id_methods)

    % Identification
    [params_id, results] = conform(sys,params_id_init,options, id_methods(i_id));   
    sys_upd = results.sys;

    % check containment of the test cases
    X0 = params_id.R0;
    for m = 1 : n_m % for each nominal trajectory

        % compute the reachable set (sometimes not computable due to 
        % by-zero-division in linerization error computation)
        params_id.u = params_true.testSuite{m}.u';
        params_id.R0 = X0 + params_true.testSuite{m}.initialState;
        %R_id = reach(sys_upd, params_id, options_reach); 

        % compute the partial reachable set Y_p (without considering
        % linerization error)
        p_GO = computeGO(sys_upd, params_id.R0.c, ...
            params_id.U.c + params_true.testSuite{m}.u', n_k);
        Y_p = cell(n_k,1);
        for k=1:n_k
            Y_p{k} = p_GO.y(:,k) + p_GO.C{k} * ...
                (params_id.R0-center(params_id.R0));
            for j = 1:k
                Y_p{k} = Y_p{k} + p_GO.D{k,j} * ...
                    (params_id.U-center(params_id.U));
            end
        end

        % check containment
        for k=1:n_k
            for s=1:length(params_true.testSuite{m}.y)
                % might fail due to numerical errors in linprog
                % (result might violate the constraints marginally, 
                % as observed for "blackNN" with "gen")
                assertLoop(contains(Y_p{k}, params_true.testSuite{m}.y(k,:,s)','exact',1e-5),i_id,m,k,s)
            end
        end
    end
end

clear mex

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
