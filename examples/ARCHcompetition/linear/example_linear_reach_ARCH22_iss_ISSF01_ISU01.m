function res = example_linear_reach_ARCH22_iss_ISSF01_ISU01()
% example_linear_reach_ARCH22_iss_ISSF01_ISU01 - ISS benchmark from ARCH22
%    solved with Alg. from [1] (constant input)
%
% Syntax:  
%    res = example_linear_reach_ARCH22_iss_ISSF01_ISU01()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References: 
%    [1] S. Bogomolov et al. "Decomposing Reach Set Computations with 
%       Low-dimensional Sets and High-Dimensional Matrices"

% Author:       Mark Wetzlinger
% Written:      09-June-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Parameters --------------------------------------------------------------

% initial set
R0 = interval(-0.0001*ones(270,1),0.0001*ones(270,1));
params.R0 = zonotope(R0);

% uncertain inputs
U  = interval([0;0.8;0.9],[0.1;1;1]);
params.U = zonotope(U);

% final time
params.tFinal = 20;


% Reachability Settings ---------------------------------------------------

options.taylorTerms   = 6;
options.zonotopeOrder = 30;
options.linAlg        = 'decomp';
options.partition     = [1, 135; 136, 270];
options.timeStep      = 0.01;


% Specification -----------------------------------------------------------

d = 5e-4;
hs1 = halfspace([0, 0, 1], -d);
hs2 = halfspace([0, 0, -1], -d);

spec = specification({hs1,hs2},'unsafeSet');


% System Dynamics ---------------------------------------------------------

% load system matrices
load iss.mat A B C

% construct the linear system object
sys = linearSys('iss',A,B,[],C);



% Reachability Analysis ---------------------------------------------------

timer = tic;
[R,res] = reach(sys,params,options,spec);
tComp = toc(timer);

disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tComp)]);



% Simulation --------------------------------------------------------------

% settings for random simulation
simOpt.points = 10;
simOpt.fracVert = 0.8;
simOpt.fracInpVert = 1;

% simulate the system
simRes = simulateRandom(sys, params, simOpt);



% Visualization -----------------------------------------------------------

figure; hold on; box on; grid on;

% plot reachable set
plotOverTime(R,3);

% plot simulation
plotOverTime(simRes,3,'Traj','y');

% plot unsafe set
plot([0 20],[d d],'--r');
plot([0 20],[-d -d],'--r');

% formatting
xlabel('t');
ylabel('y_3');
axis([0, 20, -8e-4, 8e-4]);

%------------- END OF CODE --------------