function res = test_linearSysDT_conform_07_outlier
% test_linearSysDT_conform_07_outlier - simple unit test for outlier
%   detection in reachset-conformant identification
%
% Syntax:
%    test_linearSysDT_conform_07_outlier
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Laura Luetzow, Michael Eichelbeck, Mykel Kochenderfer, and
%        Matthias Althoff. "Zono-Conformal Prediction: Zonotope-Based
%        Uncertainty Quantification for Regression and Classification
%        Tasks," arXiv, 2025.

% Authors:       Laura Luetzow
% Written:       08-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% User Specifications ----------------------------------------------------

dynamics = "pedestrian";

% specifications
n_s = 1; % not implemented for n_s > 1
n_k = 1; % not implemented for n_k > 1
n_m = 40;

% load True System Dynamics
[sys, ~,~] = loadDynamics(dynamics);
dim_u = sys.nrOfInputs;
dim_x = sys.nrOfDims;

% reachability options
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 5;
options_reach.tensorOrderOutput = 2;

% testSuite options
options_testS.p_extr = 0.1;
options_testS.stateSet = 10*interval(-ones(dim_x,1), ones(dim_x,1));

% identification options
options = options_reach;
options.cs.robustness = 1e-9;
options.cs.verbose = false;
options.cs.cost = "interval";
options.cs.constraints = "half";

% Create uncertainty set templates
c_u = zeros(dim_u,1);
c_x = zeros(dim_x,1);
G_x = [eye(dim_x) ones(dim_x,1)];
G_u = [eye(dim_u) ones(dim_u,1)];

% Initial estimates of the disturbance sets
params_true.tFinal = sys.dt * n_k - sys.dt;
params_id_init = params_true;
params_id_init.R0 = zonotope([c_x G_x]);
params_id_init.U = zonotope([c_u G_u]);

%% Create test cases
% True uncertainty sets
params_true.R0 = zonotope([0.1*randn(dim_x,1),eye(dim_x)*diag(0.01*randn(dim_u,1))]);
params_true.U = zonotope([0.1*randn(dim_u,1),eye(dim_u)*diag(rand(dim_u,1)),randn(dim_u,1)]);

testSuite = createTestSuite(sys,params_true,n_k,n_m,n_s,options_testS);

%% Identification
res = true;
params_id_init.testSuite = testSuite;

% identification results without outlier removal
[~, results_] = conform(sys,params_id_init,options);
fval0 = results_.fval;

% run all identification methods
for n_out = [1 3]
    options.cs.numOutlier = n_out;

    % test outlier detection methods "RMSE" and "MILP"
    % RMSE: remove measurements with biggest RMS error
    options.cs.outMethod = "RMSE";
    options.cs.constraints = "gen";
    [params_Rgen, results_Rgen] = conform(sys,params_id_init,options);
    numCont = aux_evaluateCoverage(sys,params_Rgen,n_k);
    assert(results_Rgen.fval <= fval0 + 1e-6) % cost should not increase
    assert(numCont >= (n_m*n_s - n_out)) % not more than n_out outliers

    % MILP: mixed-integer linear programming for optimal outlier removal
    options.cs.outMethod = "MILP";
    options.cs.constraints = "gen";
    [params_Mgen, results_Mgen] = conform(sys,params_id_init,options);
    fval_opt = results_Mgen.fval;
    numCont = aux_evaluateCoverage(sys,params_Mgen,n_k);
    assert(fval_opt <= results_Rgen.fval + 1e-6) % cost should be optimal
    assert(numCont >= (n_m*n_s - n_out)) % not more than n_out outliers

    % search and searchG: not implemented for linearSysDT
end
end


% Auxiliary functions -----------------------------------------------------

function numCont = aux_evaluateCoverage(sys,params,n_k)
X0 = params.R0;
params.tFinal = sys.dt * n_k - sys.dt;
testSuite = params.testSuite;
params = rmfield(params, 'testSuite');
numCont = 0;

for m = 1 : length(testSuite) % for each nominal trajectory

    % compute the reachable set
    params.u = testSuite(m).u;
    params.R0 = X0 + testSuite(m).x(:,1,1);

    % compute the partial reachable set Y_p
    p_GO = computeGO(sys, params.R0.c, params.U.c + testSuite(m).u, n_k);
    Y_p = cell(n_k,1);
    for k=1:n_k
        Y_p{k} = p_GO.y(:,k) + p_GO.C{k} * (params.R0-center(params.R0));
        for j = 1:k
            Y_p{k} = Y_p{k} + p_GO.D{k,j} * (params.U-center(params.U));
        end
    end

    % check containment
    for k=1:n_k
        if contains(Y_p{k}, squeeze(testSuite(m).y(:,k,:)),'exact',1e-5)
            numCont = numCont + 1;
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
