function res = example_linear_reach_ARCH22_iss_ISSC01_ISU02()
% example_linear_reach_ARCH22_iss_ISSC01_ISU02 - ISS benchmark from ARCH22
%    solved with Alg. from [1] (constant input)
%
% Syntax:  
%    res = example_linear_reach_ARCH22_iss_ISSC01_ISU02()
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

R0 = [interval(-0.0001*ones(270,1),0.0001*ones(270,1)); ...
      interval([0;0.8;0.9],[0.1;1;1])];
  
params.R0 = zonotope(R0);
params.tFinal = 20;


% Reachability Settings ---------------------------------------------------

options.taylorTerms     = 10;
options.zonotopeOrder   = 10;
options.linAlg          = 'decomp';
options.partition       = [1, 135; 136, 270; 271, 273];
options.timeStep        = 0.02;


% Spefication -------------------------------------------------------------

hs1 = halfspace([0 0 1],-1.7e-4);
hs2 = halfspace([0 0 -1],-1.7e-4);

spec = specification({hs1,hs2},'unsafeSet');


% System Dynamics ---------------------------------------------------------

% load system matrices
load iss.mat A B C

% construct extended system matrices (inputs as additional states)
dim = length(A);
A_  = [A,B;zeros(size(B,2),dim + size(B,2))];
B_  = zeros(dim+size(B,2),1);
C_  = [C,zeros(size(C,1),size(B,2))];

% construct the linear system object
sys = linearSys('iss',A_,B_,[],C_);



% Reachability Analysis ---------------------------------------------------

timer = tic;
[R,res]  = reach(sys, params, options, spec);
tComp = toc(timer);

disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tComp)]);



% Simulation --------------------------------------------------------------

% simulate the system
simRes = simulateRandom(sys, params);



% Visualization -----------------------------------------------------------

figure; hold on; 
useCORAcolors("CORA:contDynamics")

% plot reachable set
plotOverTime(R,3);

% plot initial set
plotOverTime(R(1).R0,3);

% plot simulation
plotOverTime(simRes,3);

% plot unsafe set
d1 = -hs1.d;
d2 = hs2.d;
specs = specification({
    interval([0;d1],[20;d1])
    interval([0;d2],[20;d2])
});
plot(specs, [1 2], '--')

% formatting
box on
xlabel('t');
ylabel('y_3');
axis([0, 20, -2e-4, 2e-4])
grid on

%------------- END OF CODE --------------