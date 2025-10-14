function res = test_nonlinearARX_conform_02_costs
% test_nonlinearARX_conform_02_costs - unit test for comparing
%    reachset-conformant identification with different costs
%
% Syntax:
%    res = test_nonlinearARX_conform_02_costs
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

cost_norm = {'interval', 'frob'}; 
constraints = 'half';
n_m = 2;
n_s = 100;
n_k = 4;

% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 5;
options_reach.tensorOrderOutput = 2;

% Load System Dynamics 
f = @(y,u) [0.5*y(1,1) + u(1,1) - cos(u(2,1))- 0.4*y(3,1) + u(2,1)*cos(y(1,1)); 
                0.6*y(4,1) + u(4,1)*sin(y(1,1))];
dt = 0.1;
dim_y = 2;
dim_u = 2;
p = 2;
sys = nonlinearARX('NARX2',f,dt,dim_y, dim_u, p);

params.tStart = 0;
params.tFinal = dt * (n_k-1);

% initilization
Y0{1} = zonotope([[2; 4],0.2*eye(2)]);
Y0{2} = zonotope([[2; 10],0.3*eye(2)]);
params.R0 = cartProd(Y0{1}, Y0{2});

% input
params.U = zonotope([0;0.05],0.01*eye(2));
params.u = 0.1*rand(2,n_k+1);

% Create identification data 
params.testSuite = createTestSuite(sys, params, n_k, n_m, n_s);

% Conformance Identification ----------------------------------------------

% General Options
options = options_reach;
options.cs.robustness = 1e-9;
options.cs.verbose = false;

num_id = length(cost_norm);
params_id = cell(num_id,1);
name_Rid = cell(num_id,2);

for i_id = 1:num_id
    options.cs.cost = cost_norm{i_id};
    options.cs.constraints = constraints;

    % Initial Estimates of the Disturbance Sets
    c_R0 = zeros(size(center(params.R0)));
    c_U = zeros(size(center(params.U)));
    params_id_init = params;
    params_id_init.R0 = zonotope([c_R0 eye(length(c_R0)) ones(length(c_R0),1)]);
    params_id_init.U = zonotope([c_U eye(length(c_U)) ones(length(c_U),1)]);

    % Identification ------------------------------------------------------
    [params_id{i_id}, ~] = conform(sys,params_id_init,options);   
end

% Compute Reachable Sets and Check Containment of the Test Cases
for i_id = 1 : length(params_id) % for each identification
    X0 = params_id{i_id}.R0;
    params_id{i_id}.tFinal = sys.dt * n_k - sys.dt;

    for m = 1 : n_m % for each nominal trajectory

        % compute the reachable set
        params_id{i_id}.u = params.testSuite(m).u;
        params_id{i_id}.R0 = X0 + params.testSuite(m).x(:,1,1);
        R_id = reach(sys, params_id{i_id}, options_reach);

        % compute the partial reachable set Y_p
        p_GO = computeGO(sys, params_id{i_id}.R0.c, ...
            params_id{i_id}.U.c + params.testSuite(m).u, n_k);
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
            assertLoop(contains(R_id.timePoint.set{k}, squeeze(params.testSuite(m).y(:,k,:)),'exact',1e-5),i_id,m,k)
            assertLoop(contains(Y_p{k}, squeeze(params.testSuite(m).y(:,k,:)),'exact',1e-5),i_id,m,k)
        end
    end
end

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
