function completed = example_nonlinearSysDT_conform_05_recForg
% example_nonlinearSysDT_conform_05_recForg - example for recursive 
%   reachset-conformant identification with forgetting
%
% Syntax:
%    example_nonlinearSysDT_conform_05_recForg
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
n_m = 1000;
forgetting = 0.99;
forgettingCost = 0.9; 
dynamics = "lorenz"; %"pedestrian"; 

% recursive identification methods
recMethods = ["true","whiteSlidingBS50", "recAlphaFF97","recAlphaFF90",...
    "recUnSiCoFF97","recUnSiCoFF90"];

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

% create and plot time-variant uncertainties
p_u1_all = 0.1*rand(dim_u,1)+rand(dim_u,1).*[
    zeros(dim_u,1/5*n_m) ones(dim_u,1/5*n_m) zeros(dim_u,3/5*n_m)];

figure;
sgtitle("Recursive Identification with Forgetting")
subplot(2,1,1); hold on; grid on
xlabel("$n_m$",'Interpreter','latex', 'FontSize',12)
ylabel("Uncertainties $p_{u1}$",'Interpreter','latex', 'FontSize',12)
for i=1:dim_u
    plot(1:n_m, p_u1_all(i,:));
end

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
testSuite = []; 
params_true.R0 = zonotope([0.01*randn(dim_x,1)]);
c_u_true = 0.01*randn(dim_u,1);
for m=1:n_m
    p_u1 = p_u1_all(:,m);
    params_true.U = zonotope([c_u_true,eye(dim_u)*diag(p_u1)]);
    testSuite = [testSuite; createTestSuite(sys, params_true, n_k, 1, ...
        n_s, options_testS)];
end
params_id_init.testSuite = testSuite(1:n_m);

%% Compute cost factor (necessary for computation of identification cost)
G_x_true = zeros(dim_x,0);
G_u_true = [eye(dim_u)];
cost_factors = zeros(n_m,size(G_u,2)+size(G_x,2));
cost_factors_true = zeros(n_m,size(G_u_true,2)+size(G_x_true,2));
if isa(sys,'nonlinearSysDT')
    derivatives(sys,options);
end
for m=1:n_m
    % compute the true cost for this test case
    u_nom = testSuite(m).u + c_u_true;
    x0_nom = testSuite(m).x(:,1,1) + center(params_true.R0);
    p_GO = computeGO(sys, x0_nom, u_nom, n_k);
    cost_m = 0;
    cost_m_true = 0;
    for k = 1: n_k
        % loop through time steps of the trajectory
        sum_U = 0;
        sum_U_true = 0;
        for j = 1:k
            sum_U = sum_U + abs(p_GO.D{k,j}*G_u);
            sum_U_true = sum_U_true + abs(p_GO.D{k,j}*G_u_true);
        end
        cost_m = cost_m + [abs(p_GO.C{k}*G_x) sum_U];
        cost_m_true = cost_m_true + [abs(p_GO.C{k}*G_x_true) sum_U_true];
    end
    cost_factors(m,:) = sum(cost_m,1);
    cost_factors_true(m,:) = sum(cost_m_true,1);
end

%% Identification
T = zeros(length(recMethods),1);
ms = cell(length(recMethods),1);
costs_cum = cell(length(recMethods),1);
for i_m = 1:length(recMethods)
    method = recMethods(i_m);

    % set specifications
    options.cs.forgetting = forgetting;
    options.cs.forgettingCost = forgettingCost;
    if contains(method, "FF90")
        options.cs.forgetting = 0.90;
    elseif contains(method, "FF97")
        options.cs.forgetting = 0.97;
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
    if method == "true"
        % batch identification
        T(i_m) = 0;
        ms{i_m} = 1:n_m;
        costs_cum{i_m} = zeros(1,n_m);
        for m = 1:n_m
            costs_cum{i_m}(m) = 1/m*sum(cost_factors_true(1:m,:),1)*p_u1_all(:,m);
        end
    elseif method == "whiteSlidingBS50"
        % sliding window identification
        step_sizeSliding = 50;
        ms{i_m} = step_sizeSliding:step_sizeSliding:n_m;
        params_sliding = params_id_init;
        costs_cum{i_m} = zeros(1,length(ms{i_m}));
        Ts = tic;
        m_p = 1;
        for m = ms{i_m}
            params_sliding.testSuite = testSuite(m-step_sizeSliding+1:m);
            [~, results] = conform(sys,params_sliding,options);
            costs_cum{i_m}(m_p) = 1/m*sum(cost_factors(1:m,:),1)*results.p(1:size(cost_factors,2));
            m_p = m_p +1;
        end
        T(i_m) = toc(Ts);
    else 
        % recursive identification
        ms{i_m} = step_sizeSliding:step_sizeSliding:n_m;
        costs_cum{i_m} = zeros(1,length(ms{i_m}));
        Ts = tic;
        [~, results] = conform(sys,params_id_init,options);
        T(i_m) = toc(Ts);
        m_p = 1;
        for m = ms{i_m}
            costs_cum{i_m}(m_p) = 1/m*sum(cost_factors(1:m,:),1)*results.p(:,m_p);
            m_p = m_p +1;
        end
    end
end

%% plot the results
subplot(2,1,2); hold on 
ylabel("Identification cost", 'FontSize',12)
xlabel("$n_m$",'Interpreter','latex', 'FontSize',12)
for i_m = 1:length(recMethods)
    plot(ms{i_m},costs_cum{i_m},'x-','DisplayName',recMethods(i_m),'LineWidth',1);
end
legend

completed = true;
end


% ------------------------------ END OF CODE ------------------------------
