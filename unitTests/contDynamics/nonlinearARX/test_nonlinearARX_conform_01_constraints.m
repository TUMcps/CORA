function res = test_nonlinearARX_conform_01_constraints
% test_nonlinearARX_conform_01_constraints - unit test for comparing
%        reachset-conformant identification with different constraints
%
% Syntax:
%    res = test_nonlinearARX_conform_01_constraints
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Laura Luetzow
% Written:       15-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% set random number stream
rng('default')

cost_norm = 'interval'; 
constraints = {'half', 'gen'};
n_m = 2;
n_s = 30;
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

% Reachability Settings ---------------------------------------------------

% Create identification data 
params.testSuite = createTestSuite(sys, params, n_k, n_m, n_s);

%% Conformance Identification ---------------------------------------------

% General Options
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.verbose = false;

num_id = length(constraints);
params_id = cell(num_id,1);

for i_id = 1:num_id
    options.cs.cost = cost_norm;
    options.cs.constraints = constraints{i_id};

    % Initial Estimates of the Disturbance Sets
    c_R0 = zeros(size(center(params.R0)));
    c_U = zeros(size(center(params.U)));
    params_id_init = params;
    params_id_init.R0 = zonotope([c_R0 eye(length(c_R0)) ones(length(c_R0),1)]);
    params_id_init.U = zonotope([c_U eye(length(c_U)) ones(length(c_U),1)]);

    % Identification ------------------------------------------------------
    [params_id{i_id}, ~] = conform(sys,params_id_init,options);   
end

G_X01 = params_id{1}.R0.G;
G_X02 = params_id{2}.R0.G;
G_U1 = params_id{1}.U.G;
G_U2 = params_id{2}.U.G;

tolX = max(1e-5*max([G_X01 G_X02], [],'all'), 1e-6); 
tolU = max(1e-5*max([G_U1 G_U2], [],'all'), 1e-6); 
assert(all(G_X01 <= G_X02 + tolX, 'all'))
assert(all(G_X01 >= G_X02 - tolX, 'all'))
assert(all(G_U1 <= G_U2 + tolU, 'all'))
assert(all(G_U1 >= G_U2 - tolU, 'all'))

res = true;

end


% ------------------------------ END OF CODE ------------------------------
