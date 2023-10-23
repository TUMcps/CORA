function res = example_nonlinearDA_reach_07_3bus_adaptive()
% example_nonlinearDA_reach_07_3bus_adaptive - example of nonlinear
%    differential-algebraic reachability analysis
%
% Syntax:
%    example_nonlinearDA_reach_07_3bus_adaptive()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       30-August-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Model Parameters --------------------------------------------------------

nrOfConstr = 6;
params.tFinal = 5; 

x0 = [380; 0.7];
params.y0guess = [ones(0.5*nrOfConstr, 1); zeros(0.5*nrOfConstr, 1)];
params.R0 = zonotope([x0,diag([0.1, 0.01])]);

params.U = zonotope([[1; 0.4],diag([0, 0.04])]);


% Reachability Settings ---------------------------------------------------

options.alg = 'lin-adaptive';
options.verbose = true;
options.tensorOrder = 2;


% System Dynamics ---------------------------------------------------------

sys = nonlinDASys(@bus3Dyn,@bus3Con,dim(params.R0),dim(params.U),nrOfConstr);


% Reachability Analysis ---------------------------------------------------

R = reach(sys, params, options);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.nrConstInp = 6;

simRes = simulateRandom(sys, params, simOpt);


% Visualization -----------------------------------------------------------

dim_x = [1 2];
figure; hold on; box on;

% plot reachable sets
useCORAcolors("CORA:contDynamics")
plot(R,dim_x,'DisplayName','Reachable set');

% plot initial set
plot(R(1).R0,dim_x,'DisplayName','Initial set');

% plot simulation results
plot(simRes,dim_x,'DisplayName','Simulations');

% label plot
xlabel(['x_{',num2str(dim_x(1)),'}']);
ylabel(['x_{',num2str(dim_x(2)),'}']);
legend()

% examples completed
res = true;

% ------------------------------ END OF CODE ------------------------------
