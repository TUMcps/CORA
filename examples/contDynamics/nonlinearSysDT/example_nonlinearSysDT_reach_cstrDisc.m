function completed = example_nonlinearSysDT_reach_cstrDisc
% example_nonlinearSysDT_reach_cstrDisc - example of nonlinear discrete
%    time reachability analysis, can be found in [1, Sec. 6].
% 
%
% Syntax:
%    example_nonlinearSysDT_reach_cstrDisc
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] J.M. Bravo, Robust MPC of constrained discrete-time
%        nonlinear systems based on approximated reachable sets, 2006

% Authors:       Niklas Kochdumper, Matthias Althoff
% Written:       30-January-2018
% Last update:   20-March-2020 (MA, simulateRandomDT from inherited class)
%                23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 0.15;
params.R0 = zonotope([[-0.15;-45],diag([0.005;3])]);
params.U = zonotope([zeros(2,1),diag([0.1;2])]);


% Reachability Settings  --------------------------------------------------

options.zonotopeOrder = 100;
options.tensorOrder = 3;
options.errorOrder = 5;


% System Dynamics  --------------------------------------------------------

% sampling time
dt = 0.015;

fun = @(x,u) cstrDiscr(x,u,dt);

sysDisc = nonlinearSysDT('stirredTankReactor',fun,0.015);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sysDisc,params,options);
tComp = toc;
disp("Computation time: " + tComp);


% Simulation --------------------------------------------------------------

simOpt.points = 100;
simRes = simulateRandom(sysDisc, params, simOpt);


% Visualization -----------------------------------------------------------

figure; hold on; box on;
useCORAcolors("CORA:contDynamics")

% plot reachable set
plot(R,[1 2],'DisplayName','Reachable set');

% plot initial set
plot(R.R0,[1,2],'DisplayName','Initial set');

% plot simulation
plot(simRes,[1,2],'Marker','.','LineStyle','none','DisplayName','Simulations');

% formatting
xlabel('T-T_0');
ylabel('C-C_0');
legend()

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
