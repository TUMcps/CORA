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
algorithm = {'gp','cgp'};   
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


% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 5;
options_reach.tensorOrderOutput = 2;

% Conformance Options
options = options_reach;
options.cs.robustness = 1e-9;
options.cs.verbose = false;
options.cs.cost = cost_norm;
options.cs.constraints = constraints;

% Black-box approximation options
options.id.gp_parallel = true;
options.id.gp_pop_size = 50;
options.id.gp_num_gen = 20;
options.id.cgp_num_gen = 2;
options.id.gp_parallel = false;
options.id.cgp_pop_size_base = 2;
options.id.gp_func_names = {'times','plus', 'square'};
options.id.gp_max_genes = 2;
options.id.gp_max_depth = 2;
options.id.save_res = false;
options.id.p = sys.n_p;

% Create identification data 
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s);
options.id.testSuite_id = createTestSuite(sys, params_true, ...
    n_k_train, n_m_train, n_s_train);
options.id.testSuite_val = createTestSuite(sys, params_true, ...
    n_k_val, n_m_val, n_s_val);

%% Conformance Identification ---------------------------------------------

% Initial Estimates of the Disturbance Sets
c_R0 = zeros(size(center(params_true.R0)));
c_U = zeros(size(center(params_true.U)));
params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0 eye(sys.n_p*sys.nrOfOutputs) ones(sys.n_p*sys.nrOfOutputs,1)]);
params_id_init.U = zonotope([c_U eye(sys.nrOfInputs) ...
    ones(sys.nrOfInputs,1)]);

for i_id = 1:length(algorithm)

    % Identification
    options.idAlg = algorithm{i_id};
    [params_id, results] = conform(sys,params_id_init,options, 'black');   
    sys_upd = results.sys;

    % check containment of the test cases
    X0 = params_id.R0;
    for m = 1 : n_m % for each nominal trajectory

        % compute the reachable set (sometimes not computable due to 
        % by-zero-division in linerization error computation)
        params_id.u = params_true.testSuite(m).u;
        params_id.R0 = X0 + params_true.testSuite(m).x(:,1,1);
        %R_id = reach(sys_upd, params_id, options_reach); 

        % compute the partial reachable set Y_p (without considering
        % linerization error)
        p_GO = computeGO(sys_upd, params_id.R0.c, ...
            params_id.U.c + params_true.testSuite(m).u, n_k);
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
            for s=1:length(params_true.testSuite(m).y)
                % might fail due to numerical errors in linprog
                % (result might violate the constraints marginally, 
                % as observed for "blackNN" with "gen")
                assertLoop(contains(Y_p{k}, params_true.testSuite(m).y(:,k,s),'exact',1e-5),i_id,m,k,s)
            end
        end
    end
end

clear mex

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
