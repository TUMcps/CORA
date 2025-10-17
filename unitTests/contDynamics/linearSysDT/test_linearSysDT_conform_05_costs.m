function res = test_linearSysDT_conform_05_costs
% test_linearSysDT_conform_05_costs - unit test for comparing
%        reachset-conformant identification with different costs
%
% Syntax:
%    res = test_linearSysDT_conform_05_costs
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       08-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set random number stream
rng('default')

cost_norm = {"interval","interval1","interval5","frob"}; 
n_m = 2;
n_s = 5;
n_k = 3;

% Reachability Settings  
options_reach.zonotopeOrder = 100;

% Load System Dynamics 
[sys, params_true.R0, params_true.U] = loadDynamics('pedestrian');
params_true.tFinal = sys.dt * n_k - sys.dt;

% Create identification data 
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s);

% Conformance Identification ----------------------------------------------

% General Options
options = options_reach;
options.cs.robustness = 1e-9;
options.cs.verbose = false;


for constraints = ["gen","half"]
    options.cs.constraints = constraints;
    num_id = length(cost_norm);
    params_id = cell(num_id,1);
    for i_id = 1:num_id
        if strcmp(cost_norm{i_id},'interval1')
            % one rotation
            options.cs.cost = 'interval';
            options.cs.numRotations = 1;
        elseif strcmp(cost_norm{i_id},'interval5')
            % five rotations
            options.cs.cost = 'interval';
            options.cs.numRotations = 5;
        else
            % no rotations
            options.cs.cost = cost_norm{i_id};
            options.cs.numRotations = 0;
        end

        % Initial Estimates of the Disturbance Sets
        c_R0 = zeros(size(center(params_true.R0)));
        c_U = zeros(size(center(params_true.U)));
        params_id_init = params_true;
        params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfDims) ones(sys.nrOfDims,1)]);
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
            params_id{i_id}.u = params_true.testSuite(m).u;
            params_id{i_id}.R0 = X0 + params_true.testSuite(m).x(:,1,1);
            R_id = reach(sys, params_id{i_id}, options_reach);

            % compute the partial reachable set Y_p
            p_GO = computeGO(sys, params_id{i_id}.R0.c, ...
                params_id{i_id}.U.c + params_true.testSuite(m).u, n_k);
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
                assertLoop(contains(R_id.timePoint.set{k}, squeeze(params_true.testSuite(m).y(:,k,:)),'exact',1e-5),i_id,m,k)
                assertLoop(contains(Y_p{k}, squeeze(params_true.testSuite(m).y(:,k,:)),'exact',1e-5),i_id,m,k)
            end
        end
    end
end

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
