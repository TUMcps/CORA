function res = testLong_nonlinearSysDT_conform_03_gray
% testLong_nonlinearSysDT_conform_03_gray - unit test for computing 
%   rechset-conformant scaling factor and center vectors
%
% Syntax:
%    res = testLong_nonlinearSysDT_conform_03_gray
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

dyn = 'lorenz_2D';
cost_norm = "interval"; 
constraints = "half";
methodsGray = ["graySeq","grayLS","graySim"];
n_m = 2;
n_s = 30;
n_k = 3;

% Reachability Settings  
options_reach.verbose = false;
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 5;
options_reach.tensorOrderOutput = 2;

% Load System Dynamics 
[sys, params_true.R0, params_true.U, p_true] = loadDynamics(dyn,'diag');
params_true.tFinal = sys.dt * n_k - sys.dt;

% Create identification data 
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s);

%% Conformance Identification ---------------------------------------------


% General Options
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.verbose = false;
options.cs.cost = cost_norm;
options.cs.constraints = constraints;

% Initial Estimates of the Disturbance Sets
c_R0 = zeros(size(params_true.R0.c));
c_U = zeros(size(params_true.U.c));
params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfDims)]);
params_id_init.U = zonotope([c_U eye(sys.nrOfInputs)]);

for val = 1:2
    % first iteration: use default set_p function (= only estimate the
    % center vectors)

    % second iteration: use user-defined set_p function (= estimate some
    % parameters and the center vector of U)
    if val == 2
        options.cs.p0 = 0.01*randn(size(p_true,1)+sys.nrOfInputs,1); % estimate parameter and center of U
        options.cs.set_p = @(p,params) aux_set_p(p, params, dyn);
        options.cs.derivRecomputation = false;
    end

    for type = methodsGray
        % Identification ------------------------------------------------------
        [params_id, results] = conform(sys,params_id_init,options,type);
        sys_id = results.sys;

        X0 = params_id.R0;
        params_id.tFinal = sys.dt * n_k - sys.dt;
        params_id = rmfield(params_id, 'testSuite');

        for m = 1 : n_m % for each nominal trajectory

            % compute the reachable set
            params_id.u = params_true.testSuite{m}.u';
            params_id.R0 = X0 + params_true.testSuite{m}.initialState;
            R_id = reach(sys_id, params_id, options_reach);

            % compute the partial reachable set Y_p
            p_GO = computeGO(sys_id, params_id.R0.c, ...
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
                    assertLoop(contains(R_id.timePoint.set{k}, params_true.testSuite{m}.y(k,:,s)','exact',1e-5),val,type,m,k,s)
                    assertLoop(contains(Y_p{k}, params_true.testSuite{m}.y(k,:,s)','exact',1e-5),val,type,m,k,s)
                end
            end
        end
    end
end
res = true;
end


% Auxiliary functions -----------------------------------------------------

function [sys,params] = aux_set_p(p, params, dyn)
% return dynamics and uncertainty sets parameterized by p
[sys, ~, ~, p_t] = loadDynamics(dyn, 'diag', p);
c = p(length(p_t)+1:end);
params.U = zonotope(c, params.U.generators);
end

% ------------------------------ END OF CODE ------------------------------
