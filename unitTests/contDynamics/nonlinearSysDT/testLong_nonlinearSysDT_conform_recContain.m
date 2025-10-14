function res = testLong_nonlinearSysDT_conform_recContain
% testLong_nonlinearSysDT_conform_recContain - unit test for recursive
%   reachset-conformant identification with containment check
%
% Syntax:
%    testLong_nonlinearSysDT_conform_recContain
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1]

% Authors:       Laura Luetzow
% Written:       04-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% User Specifications ----------------------------------------------------

dynamics = "lorenz";

% specifications
n_s = 1;
n_k = 4;
n_m = 70;
numPoints = 1;
numCluster = 10;
batch_size = 20;
forgetting = 1;
forgettingCost = 1;

% recursive identification methods
recMethods = ["white", "recAlpha", "recAlphaBS10", "recAlphaFF90", ...
    "recUnSiCo", "recUnSiCoQ5", "recUnSiCoFF90"];

% Set random number stream
rng(1)

% load True System Dynamics
[sys, ~,~] = loadDynamics(dynamics);
dim_u = sys.nrOfInputs;
dim_x = sys.dim;

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

% run all identification methods
for i_id = 1:length(recMethods)
    method = recMethods(i_id);

    % set specifications
    options.cs.batchSize = batch_size;
    options.cs.numPoints = numPoints;
    options.cs.numCluster = numCluster;
    options.cs.forgetting = forgetting;
    options.cs.forgettingCost = forgettingCost;

    % index of first test case that must be reachset-conformant
    m_0 = 1;

    % set recursive identification method
    if contains(method, "recUnSiCo")
        options.cs.recMethod = "recUnSiCo";
    elseif contains(method, "recAlpha")
        options.cs.recMethod = "recAlpha";
    else
        options.cs.recMethod = "";
    end

    % specify number of points for UnSiCo
    if contains(method, "Q5")
        options.cs.numPoints = 5;
    end

    % specify batch size
    if contains(method, "BS10")
        options.cs.batchSize = 10;
    end

    % specify forgetting factor
    if contains(method, "FF90")
        options.cs.forgetting = 0.90;
        options.cs.forgettingCost = 0.90;
        num_recursions = ceil(n_m/options.cs.batchSize);
        m_0 = (num_recursions-1)*options.cs.batchSize + 1;
    end

    % identification
    [params_id, results] = conform(sys,params_id_init,options);
    if isempty(results.fval) || isempty(results.p) || any(~isfinite(results.fval(end)),'all') || any(~isfinite(results.p),'all')
        res = false;
    end

    % check containment
    X0 = params_id.R0;
    params_id.tFinal = sys.dt * n_k - sys.dt;
    params_id = rmfield(params_id, 'testSuite');

    for m = m_0 : n_m % for each nominal trajectory

        % compute the reachable set
        params_id.u = params_id_init.testSuite(m).u;
        params_id.R0 = X0 + params_id_init.testSuite(m).x(:,1,1);
        R_id = reach(sys, params_id, options_reach);

        % check containment
        for k=1:n_k
            assertLoop(contains(R_id.timePoint.set{k}, squeeze(params_id_init.testSuite(m).y(:,k,:)),'exact',1e-4),i_id,m,k)
        end
    end
end
end


% ------------------------------ END OF CODE ------------------------------
