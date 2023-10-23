function completed = example_nonlinearDA_reach_01_powerSystem_3bus()
% example_nonlinearDA_reach_01_powerSystem_3bus - example of 
%    nonlinear-differential-algebraic reachability analysis, 3-bus power system
%
% Syntax:
%    completed = example_nonlinearDA_reach_01_powerSystem_3bus()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false

% Authors:       Matthias Althoff
% Written:       18-August-2016
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameter ---------------------------------------------------------------

nrOfConstr = 6;
params.tFinal = 5; 

x0 = [380; 0.7];
params.y0guess = [ones(0.5*nrOfConstr, 1); zeros(0.5*nrOfConstr, 1)];
params.R0 = zonotope([x0,diag([0.1, 0.01])]);

params.U = zonotope([[1; 0.4],diag([0, 0.04])]);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.05;
options.taylorTerms = 6;
options.zonotopeOrder = 200;
options.tensorOrder = 2;
options.errorOrder = 10;

options.maxError_x = [0.5; 0];
options.maxError_y = 0.005*[1; 1; 1; 1; 1; 1];


% System Dynamics ---------------------------------------------------------

powerDyn = nonlinDASys(@bus3Dyn,@bus3Con);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(powerDyn, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 10;
simRes = simulateRandom(powerDyn, params, simOpt);


% Visualization -----------------------------------------------------------

projDim = [1 2];
    
figure; hold on; box on;

% plot reachable sets
useCORAcolors("CORA:contDynamics")
plot(R,projDim,'DisplayName','Reachable set');

% plot initial set
plot(R.R0,projDim, 'DisplayName','Initial set');

% plot simulation results      
plot(simRes,projDim,'DisplayName','Simulations');

% label plot
xlabel(['x_{',num2str(projDim(1)),'}']);
ylabel(['x_{',num2str(projDim(2)),'}']);

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
