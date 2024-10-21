function res = test_nonlinearSysDT_conform_02_costs
% test_nonlinearSysDT_conform_02_costs - unit test for comparing
%        reachset-conformant identification with different costs
%
% Syntax:
%    res = test_nonlinearSysDT_conform_02_costs
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

cost_norm = {"interval", "frob"}; 
constraints = "half";
n_m = 2;
n_s = 100;
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

% Conformance Identification ----------------------------------------------

% General Options
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.verbose = false;

num_id = length(cost_norm);
params_id = cell(num_id,1);
options.cs.constraints = constraints;

for i_id = 1:num_id
    options.cs.cost = cost_norm{i_id};

    % Initial Estimates of the Disturbance Sets
    c_R0 = zeros(size(center(params_true.R0)));
    c_U = zeros(size(center(params_true.U)));
    params_id_init = params_true;
    params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfStates) ones(sys.nrOfStates,1)]);
    params_id_init.U = zonotope([c_U eye(sys.nrOfInputs) ...
        ones(sys.nrOfInputs,1)]);

    % Identification ------------------------------------------------------
    [params_id{i_id}, ~] = conform(sys,params_id_init,options);   
end

% Compute Reachable Sets and Check Containment of the Test Cases
for i_id = 1 : length(params_id) % for each identification
    X0 = params_id{i_id}.R0;
    params_id{i_id}.tFinal = sys.dt * n_k - sys.dt;
    params_id{i_id} = rmfield(params_id{i_id}, 'testSuite');

    for m = 1 : n_m % for each nominal trajectory

        % compute the reachable set
        params_id{i_id}.u = params_true.testSuite{m}.u';
        params_id{i_id}.R0 = X0 + params_true.testSuite{m}.initialState;
        R_id = reach(sys, params_id{i_id}, options_reach);

        % compute the partial reachable set Y_p
        p_GO = computeGO(sys, params_id{i_id}.R0.c, ...
            params_id{i_id}.U.c + params_true.testSuite{m}.u', n_k);
        Y_p = cell(n_k,1);
        for k=1:n_k
            Y_p{k} = p_GO.y(:,k) + p_GO.C{k} * ...
                (params_id{i_id}.R0-center(params_id{i_id}.R0));
            for j = 1:k
                Y_p{k} = Y_p{k} + p_GO.D{k,j} * ...
                    (params_id{i_id}.U-center(params_id{i_id}.U));
            end
        end

        % check containment
        for k=1:n_k
            [~,n,numPoints] = size(params_true.testSuite{m}.y);
            assertLoop(contains(R_id.timePoint.set{k}, reshape(params_true.testSuite{m}.y(k,:,:),n,numPoints,1),'exact',1e-5),i_id,m,k)
            assertLoop(contains(Y_p{k}, reshape(params_true.testSuite{m}.y(k,:,:),n,numPoints,1),'exact',1e-5),i_id,m,k)
        end
    end
end

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
