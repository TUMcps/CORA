function res = example_linear_reach_ARCH21_beam_CB01()
% example_linear_reach_ARCH21_beam_CB01 - clamped beam benchmark
%    from the 2021 ARCH competition
%
% Syntax:  
%    res = example_linear_reach_ARCH21_beam_CB01()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 

% Author:       Mark Wetzlinger
% Written:      02-June-2021
% Last update:  17-June-2021 (damping case)
% Last revision:---

%------------- BEGIN CODE --------------

sys_conv = beam();

res = 1;

% Model Derivation --------------------------------------------------------

% nodes in model
N = 100;
% node of interest
node = 70;

% constants
rho = 7.3e-4;   % density
L = 200;        % length of beam
Q = 1;          % cross-section area (renamed from A)
E = 30e6;       % Young's modulus

ell = L/N;      % length of individual discrete element

% mass matrix (NxN)
M = (rho*Q*ell) / 2 * diag([2*ones(N-1,1);1]);
Minv = M^(-1);

% load
F = zonotope(10000,100);

% tridiagonal matrix (NxN)
mat = zeros(N);
mat(1,1) = 2; mat(1,2) = -1;
mat(N,N-1) = -1; mat(N,N) = 1;
for r=2:N-1
    mat(r,1+(r-2)) = -1;
    mat(r,2+(r-2)) = 2;
    mat(r,3+(r-2)) = -1;
end
% stiffness matrix (NxN)
K = E*Q/ell * mat;
% damping matrix (NxN)
a = 1e-6;
b = 1e-6;
D = a*K + b*M;

% state matrix (no damping)
% A = [zeros(N) eye(N); -Minv*K zeros(N)];
% state matrix (damped)
A = [zeros(N) eye(N); -Minv*K -Minv*D];

% Parameters --------------------------------------------------------------

params.tFinal = 0.01;

% nr of states
dim_x = length(A);

% initial set: bar at rest
X0_C = cartProd( zonotope(zeros(dim_x,1)), Minv(end,end)*F );
X0_F = zonotope(zeros(dim_x,1));

% input set
U_C = zonotope(0);
U_F = cartProd( zonotope(zeros(dim_x-1,1)), Minv(end,end)*F );


% Reachability Settings ---------------------------------------------------

options.taylorTerms = 20;
options.zonotopeOrder = 25;

linAlg_C = 'standard';
linAlg_F = 'wrapping-free';

timeStep_C = params.tFinal/7000;
timeStep_F = params.tFinal/10000;


% System Dynamics ---------------------------------------------------------

% 1. constant inputs (extended state matrix)
A_C = [ A, [zeros(dim_x-1,1);1]; zeros(1,dim_x+1) ];
B_C = 0;
% C_C = zeros(2,2*N+1);
% C_C(1,node) = 1; C_C(2,2*node) = 1;

% 2. time-varying inputs
A_F = A;
B_F = 1; % (canonical form)
% C_F = zeros(2,2*N);
% C_F(1,node) = 1; C_F(2,2*node) = 1;


% construct linear system objects
sys_C = linearSys(['CB21C_' num2str(N)],A_C,B_C);
sys_F = linearSys(['CB21F_' num2str(N)],A_F,B_F);


% Reachability Analysis ---------------------------------------------------

% 1. constant inputs
timer = tic;

params.R0 = X0_C;
params.U = U_C;
options.timeStep = timeStep_C;
options.linAlg = linAlg_C;
R_C = reach(sys_C, params, options);

tComp_C = toc(timer);
disp(['computation time (CB21C_' num2str(N) '): ' num2str(tComp_C)]);

[x_max_C,v_max_C] = maxvals(R_C,node);
disp(['x_' num2str(node) ',max (C): ' num2str(x_max_C)]);
disp(['v_' num2str(node) ',max (C): ' num2str(v_max_C)]);

% 2. time-varying inputs
timer = tic;

params.R0 = X0_F;
params.U = U_F;
options.timeStep = timeStep_F;
options.linAlg = linAlg_F;
R_F = reach(sys_F, params, options);

tComp_F = toc(timer);
disp(['computation time (CB21F_' num2str(N) '): ' num2str(tComp_F)]);

[x_max_F,v_max_F] = maxvals(R_F,node);
disp(['x_' num2str(node) ',max (F): ' num2str(x_max_F)]);
disp(['v_' num2str(node) ',max (F): ' num2str(v_max_F)]);


% Simulation --------------------------------------------------------------

% 1. constant inputs
params.R0 = X0_C;
params.U = U_C;
% settings for random simulation
simOpt.points = 2; % only two extreme trajectories since no inputs
simOpt.fracVert = 1;
simOpt.fracInpVert = 0; % no inputs
simOpt.inpChanges = 0;  % no inputs

simRes_C = simulateRandom(sys_C, params, simOpt);

% 2. time-varying inputs
params.R0 = X0_F;
params.U = U_F;
% settings for random simulation
simOpt.points = 20;
simOpt.fracVert = 0.1;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 10;

simRes_F = simulateRandom(sys_F, params, simOpt);


% Visualization -----------------------------------------------------------

% requirement: velocity at node 70/100 at t \in [8.15e-3, 8.4e-3]
projDim = node * 2;

% 1. constant inputs 
figure; hold on; box on;
title("Constant uncertain inputs");
plotOverTime(R_C,projDim,'FaceColor',[.6 .6 .6],'EdgeColor','none');
plotOverTime(simRes_C,projDim);
xlim([8.15e-3,8.4e-3]);
xlabel('t');
ylabel(['v_{' num2str(node) '}']);

% 2. time-varying inputs 
figure; hold on; box on;
title("Time-varying inputs");
plotOverTime(R_F,projDim,'FaceColor',[.6 .6 .6],'EdgeColor','none');
plotOverTime(simRes_F,projDim);
xlim([8.15e-3,8.4e-3]);
xlabel('t');
ylabel(['v_{' num2str(node) '}']);

end


% Auxiliary Function
function [x,v] = maxvals(R,N)

% init maxima
x = -Inf;
v = -Inf;

% maximum value at node N (x_N and x_2N for position and velocity resp.)
for k=1:length(R.timeInterval.set)
    temp = interval(project(R.timeInterval.set{k},[N,2*N]));
    supk = supremum(temp);
    if supk(1) > x
        x = supk(2);
    end
    if supk(2) > v
        v = supk(2);
    end
end

end

%------------- END OF CODE --------------