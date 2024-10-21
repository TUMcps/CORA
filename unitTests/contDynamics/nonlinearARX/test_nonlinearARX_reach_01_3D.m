function res = test_nonlinearARX_reach_01_3D
% test_nonlinearARX_reach_01_3D - test of nonlinear reachability 
%    analysis for NARX models.
%
% Syntax:
%    completed = test_nonlinearARX_reach_01_3D
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
    0.4*y(3,1) + u(2,1)*cos(y(1,1)); 0.6*y(5,1) + u(4,1)*sin(y(1,1))];
dt = 0.1;
dim_y = 3;
dim_u = 2;
p = 2;
sys = nonlinearARX(f,dt,dim_y, dim_u, p);

% Parameters --------------------------------------------------------------

% time horizon
N_k = 20;
params.tStart = 0;
params.tFinal = dt * (N_k-1);

% initilization
params.Y0{1} = zonotope([[2; 4; 4],0.2*eye(3)]);
params.Y0{2} = zonotope([[2; 10; 4],0.3*eye(3)]);


% input
params.U = zonotope([0;0.05],0.01*eye(2));
params.u = 0.1*rand(2,N_k+1);

% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 50; %zonotope order
options.tensorOrder = 2;
options.errorOrder = 1;
options.lagrangeRem.simplify = 'simplify';


% Reachability Analysis ---------------------------------------------------

options.armaxAlg = 'exactAddition';
R_exact = reach(sys, params, options);
options.armaxAlg = 'tvpGeneral';
R_gen = reach(sys, params, options);


% Simulation --------------------------------------------------------------

sim_points = 30;
y_init1 = randPoint(params.Y0{1},sim_points, 'extreme');
y_init2 = randPoint(params.Y0{2},sim_points, 'extreme');
u = params.u;
y_sim = cell(sim_points,1);

for i = 1:sim_points
    params.y0 = [y_init1(:,i) y_init2(:,i)];
    params.u = u + randPoint(params.U, size(params.u, 2), 'extreme');
    [~,~,~,y_sim{i}] = simulate(sys, params);
end
    
% check equality of the reachable sets and containment of the sample points
for k = 1: length(R_exact.timePoint.set)
    Yk_exact = zonotope(R_exact.timePoint.set{k});
    Yk_gen = zonotope(R_gen.timePoint.set{k});
    assertLoop(isequal(reduce(Yk_exact,'girard',2),reduce(Yk_gen,'girard',2),1e-3),k)
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
