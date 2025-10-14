function completed = example_nonlinearSysDT_conform_04_rec
% example_nonlinearSysDT_conform_04_rec - example for recursive 
%   reachset-conformant identification
%
% Syntax:
%    example_nonlinearSysDT_conform_04_rec
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] 

% Authors:       Laura Luetzow
% Written:       04-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% User Specifications ----------------------------------------------------

% specifications
n_s = 1; 
n_k = 4;
n_m_sizes = [50 [1:5]*200];
numPoints = 1;
numCluster = 10;
batch_size = 50;
dynamics = "lorenz"; %"pedestrian"; 

% recursive identification methods
recMethods = ["white", "recAlphaBS10", "recAlphaBS50", ...
    "recUnSiCoQ1", "recUnSiCoQ5"];

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
testSuite = createTestSuite(sys, params_true, n_k, max(n_m_sizes), ...
    n_s, options_testS);

%% Identification
T = zeros(length(recMethods),length(n_m_sizes));
costs_cum = zeros(length(recMethods),length(n_m_sizes));
for i_nm = 1:length(n_m_sizes)
    n_m_i = n_m_sizes(i_nm);
    params_id_init.testSuite = testSuite(1:n_m_i);

    for i_m = 1:length(recMethods)
        method = recMethods(i_m);
        
        % set specifications
        options.cs.batchSize = batch_size;
        options.cs.numPoints = numPoints;
        options.cs.numCluster = numCluster;
        if contains(method, "Q1")
            options.cs.numPoints = 1;
        elseif contains(method, "Q5")
            options.cs.numPoints = 5;
        end
        if contains(method, "BS50")
            options.cs.batchSize = 50;
        elseif contains(method, "BS10")
            options.cs.batchSize = 10;
        end

        % set recursive identification method
        if contains(method, "recUnSiCo")
            options.cs.recMethod = "recUnSiCo";
        elseif contains(method, "recAlpha")
            options.cs.recMethod = "recAlpha";
        else
            options.cs.recMethod = "";
        end

        % identification
        Ts = tic;
        [~, results] = conform(sys,params_id_init,options);
        T(i_m,i_nm) = toc(Ts);
        costs_cum(i_m,i_nm) = 1/n_m_i*results.fval(end);
    end
end

%% plot the results
figure

% plot identification cost
subplot(2,1,1); hold on 
sgtitle("Recursive Identification")
ylabel("Identification cost", 'FontSize',12)
xlabel("$n_m$",'Interpreter','latex', 'FontSize',12)
for i_m = 1:length(recMethods)
    plot(n_m_sizes,costs_cum(i_m,:),'x-','DisplayName',recMethods(i_m),'LineWidth',1);
end

% plot computation time
subplot(2,1,2); hold on 
ylabel("Computation time [s]", 'FontSize',12)
xlabel("$n_m$",'Interpreter','latex', 'FontSize',12)
for i_m = 1:length(recMethods)
    plot(n_m_sizes,T(i_m,:),'x-','DisplayName',recMethods(i_m),'LineWidth',1);
end
legend

completed = true;
end


% ------------------------------ END OF CODE ------------------------------
