function completed = example_nonlinear_reach_16_infiniteBusRoA_polynomialized
% example_nonlinear_reach_16_infiniteBusRoA_polynomialized - example for 
%    verifying the region of attraction of a single-machine-infinite-bus 
%    system from [1] after polynomialization
%
% Syntax:
%    completed = example_nonlinear_reach_16_infiniteBusRoA_polynomialized
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% References:
%    [1] M. Althoff. "Formal Verification of Power Systems", submitted to
%        ARCH 2022.

% Authors:       Matthias Althoff
% Written:       02-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 1;
params.U = zonotope([0, 0]);

% over-approximate initial set using interval arithmetic
R0_origin = zonotope([[asin(1/5); 0],0.05*eye(2)]);
R0_interval = interval(R0_origin);
R0_trans = [...
    sin(R0_interval(1));...
    1 - cos(R0_interval(1));...
    R0_interval(2)];
params.R0 = zonotope(R0_trans);


% Reachability Settings ---------------------------------------------------

% settings
options.timeStep = 0.0005;
options.taylorTerms = 4;
options.zonotopeOrder = 1000;
options.alg = 'poly';
options.tensorOrder = 3;
options.errorOrder = 2;
options.intermediateOrder = 2;

% settings for polynomial zonotopes
polyZono.maxDepGenOrder = 30;
polyZono.maxPolyZonoRatio = 0.001;
polyZono.restructureTechnique = 'reduceFullGirard';
options.polyZono = polyZono;

options.lagrangeRem.simplify = 'simplify';


% System Dynamics ---------------------------------------------------------

sys = nonlinearSys(@infiniteBus_polynomialized);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simRes = simulateRandom(sys, params, simOpt);


% Visualization -----------------------------------------------------------

figure; hold on; box on;
useCORAcolors("CORA:contDynamics")

% plot reachable sets
plot(R);

% plot initial set
plot(R.R0,[1 2]);

% plot simulation results     
plot(simRes);

% label plot
xlabel('x_{1}');
ylabel('x_{2}');


% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
