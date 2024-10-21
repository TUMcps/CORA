function res = test_linearARX_reach
% test_linearARX_reach - unit test for reach where the reachable set 
%   of ARX models is computed with different algorithmen
%
% Syntax:
%    res = test_linearARX_reach
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Example: 
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.

% Authors:       Laura Luetzow
% Written:       03-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

num_samples = 200;

% set random number stream
rng('default')

%% define linearSys -------------------------------------------------------

% stable system matrix: n x n
A = [0.9618    0.0286    0.0509   -0.0265
    0.0179    1.0257    0.0548    0.0836
    0.0066   -0.0111    0.9378    0.0365
    0.0990   -0.0478    0.0811    0.9747];
dim_x = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
    2 1 0;
    0 0 1;
    0 -2 1];
dim_u = size(B,2);

% output matrix: q x n
C = [1 1 0 0;
    0 -0.5 0.5 0];

% throughput matrix: q x m
D = [0 0 1;
    0 0 0];
dim_y = size(D,1);

% disturbance matrix: n x r
E = eye(dim_x);

% noise matrix: q x s
F = eye(dim_y);

% time step
dt = 0.1;

% initialize system objects
sys_lin = linearSysDT(A,B,[],C,D,[],E,F,dt);
sys_ARX = linearARX(sys_lin);
p = sys_ARX.n_p;

%% model parameters and reachability settings -----------------------------
% general parameters ------------------------------------------------------

% disturbance sets
U = zonotope(zeros(dim_u,1));
W = 10*zonotope(0.02+zeros(dim_x,1),0.02*diag(ones(dim_x,1)));
V = 10*zonotope(-0.01+zeros(dim_y,1),0.01*[diag(ones(dim_y,1)) ones(dim_y,1)]);

% true initial state (unknown)
x0 = 10*ones(dim_x,1);

% time horizon
dt_steps = 11;
tFinal = dt * dt_steps;
tStart = 0;

% input and distuabance vectors
u = randn(dim_u,dt_steps+1);
w = randPoint(W,dt_steps);
v = randPoint(V,dt_steps+1);

%% initialization ---------------------------------------------------------

% compute the initial measurements y(1),..., y(p) -------------------------
% paramters for computing the initial measurements
params_init.x0 = x0;
params_init.u = u(:,1:p);
params_init.w = w(:,1:p-1);
params_init.v = v(:,1:p);
params_init.tFinal = (p-1)*dt;

% compute the initial measurements
[~,~,~,y_init] = simulate(sys_lin,params_init);


%% simulation of sample points --------------------------------------------

% parameters for simulating the linearARX
params_ARX_sim.tStart = tStart;
params_ARX_sim.tFinal = tFinal;
params_ARX_sim.y0 = y_init';

% simulate initial point
y_sim = cell(num_samples,1);
setTVP(sys_ARX);

% simulate random points
for i_sample=1:num_samples
    if i_sample < num_samples/2
        % 300 points with maximum disturbances
        v = randPoint(V,dt_steps+1, 'extreme');
        w = randPoint(W,dt_steps, 'extreme');
    else 
        % random points
        v = randPoint(V,dt_steps+1);
        w = randPoint(W,dt_steps);
    end
    params_ARX_sim.u = [u; w zeros(dim_x,1); v];
    [~,~,~,y_sim{i_sample}] = simulate(sys_ARX,params_ARX_sim);
end

%% reachability analysis --------------------------------------------------

% reachability parameters for linearARX
params_ARX_reach.tStart = tStart;
params_ARX_reach.tFinal = tFinal;
params_ARX_reach.y0 = y_init';
params_ARX_reach.U = cartProd(U, cartProd(W,V));
params_ARX_reach.u = [u; zeros(dim_x+dim_y, size(u,2))];

% compute the reachable set with ARX & the algorithms from [1]
ind = 1;
for reach_alg = {'exactAddition' 'tvpGeneral' 'tvpEfficient'}
    options.armaxAlg = reach_alg{1};
    R_ARX{ind} = reach(sys_ARX,params_ARX_reach,options);
    R_ARX{ind} = R_ARX{ind}.timePoint.set;

    % check if all measurements are contained in the reachable set of the 
    % ARX model
    assert(aux_checkContainment(R_ARX{ind}, y_sim, p+1));

    % check if the results of the different ARX reachability algorithms 
    % are equal
    if ind >= 2
        assert(aux_checkEquality(R_ARX{ind-1}, R_ARX{ind}, p+1));
    end
    ind = ind + 1;
end
res = true;
end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkContainment(R, y, k_start)
% check if the sample points y are contained in the set R

tol_contains = 1e-9;
res = true;
for i_y = 1:length(y)
    y_i = y{i_y};
    for k=k_start:length(R)
        % check containment
        assertLoop(contains(R{k}, y_i(k,:)', 'exact', tol_contains),i_y,k)
    end
end
end

function res = aux_checkEquality(R1, R2, k_start)
% check equality of the reachable sets R1 and R2

tol_iseuqal = 1e-9;
res = true;

for k=k_start:length(R1)
    % check equality
    assertLoop(isequal(R1{k}, R2{k}, tol_iseuqal),k)
end
end

% ------------------------------ END OF CODE ------------------------------
