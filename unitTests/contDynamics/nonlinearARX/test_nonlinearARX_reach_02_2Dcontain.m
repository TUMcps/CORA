function res = test_nonlinearARX_reach_02_2Dcontain
% test_nonlinearARX_reach_02_2Dcontain - test of nonlinear reachability 
%    analysis for NARX models.
%
% Syntax:
%    completed = test_nonlinearARX_reach_02_2Dcontain
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 

% Authors:       Laura Luetzow
% Written:       13-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% set random number stream
rng('default')

% System Dynamics ---------------------------------------------------------

f = @(y,u) [0.5*y(1,1) + u(1,1) - cos(u(2,1)); ...
    0.4*y(3,1) + u(5,1)*cos(y(4,1)) + 0.6*y(5,1) + u(7,1)*sin(y(1,1))];
dt = 0.1;
dim_y = 2;
dim_u = 2;
p = 3;
sys = nonlinearARX(f,dt,dim_y, dim_u, p);

% Parameters --------------------------------------------------------------

% time horizon
N_k = 20;
params.tStart = 0;
params.tFinal = dt * (N_k-1);

% initilization
Y0{1} = zonotope([[2; 4],0.2*eye(2)]);
Y0{2} = zonotope([[-1; 10],diag([0.5 0.1])]);
Y0{3} = zonotope([[5; 2],diag([0.05 0.3])]);
params.R0 = getR0(sys, Y0);

% input
params.U = zonotope([0;0.05],0.01*eye(2));
params.u = 0.1*rand(2,N_k+1);

% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 50; %zonotope order
options.tensorOrder = 2;
options.errorOrder = 1;
options.lagrangeRem.simplify = 'simplify';


% Reachability Analysis ---------------------------------------------------

tic
options.armaxAlg = 'exactAddition';
R_exact = reach(sys, params, options);
options.armaxAlg = 'tvpGeneral';
R_gen = reach(sys, params, options);


% Simulation --------------------------------------------------------------

sim_points = 50;
x_init1 = randPoint(params.R0,sim_points, 'extreme');
u = params.u;
tol_contains = 1e-9;
y_sim = cell(sim_points,1);

for i = 1:sim_points
    params.x0 = x_init1(:,i);
    params.u = u + randPoint(params.U, size(params.u, 2), 'extreme');
    [tVec,~,~,y_sim{i}] = simulate(sys, params);
end
    
% check equality of the reachable sets and containment of the sample points
for k = 1: length(R_exact.timePoint.set)
    Yk_exact = zonotope(R_exact.timePoint.set{k});
    Yk_gen = zonotope(R_gen.timePoint.set{k});
    assertLoop(isequal(reduce(Yk_exact,'girard',2),reduce(Yk_gen,'girard',2),1e-3),k)

    for i = 1:sim_points
        assertLoop(contains(Yk_exact, y_sim{i}(k,:)', 'exact', tol_contains),k,i)
        assertLoop(contains(Yk_gen, y_sim{i}(k,:)', 'exact', tol_contains),k,i)
    end
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
