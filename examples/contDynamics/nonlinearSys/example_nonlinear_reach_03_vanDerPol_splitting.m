function completed = example_nonlinear_reach_03_vanDerPol_splitting()
% example_nonlinear_reach_03_vanDerPol_splitting - example of nonlinear
%    reachability analysis with splitting, can be found in [1, Sec. 3.4.5]
%    or in [2]. A new technique for computing this example with less
%    spliiting has been published in [3].
%
% Syntax:
%    completed = example_nonlinear_reach_03_vanDerPol_splitting()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010
%    [2] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008
%    [3] M. Althoff. "Reachability analysis of nonlinear systems using 
%        conservative polynomialization and non-convex sets", HSCC 2013

% Authors:       Matthias Althoff
% Written:       26-June-2009
% Last update:   23-April-2020 (restucture params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 0.5;

% specify initial set as zonotope bundle so that splitting sets is exact
Z0{1} = zonotope([1.4 0.3 0; 2.3 0 0.05]);
params.R0 = zonoBundle(Z0);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.02;
options.taylorTerms = 4;
options.zonotopeOrder = 10;
options.intermediateOrder = 10;
options.errorOrder = 5;

options.alg = 'lin';
options.tensorOrder = 3;

options.maxError = 0.05*[1; 1];
options.reductionInterval = 100;
options.verbose = true;


% System Dynamics ---------------------------------------------------------

vanderPol = nonlinearSys(@vanderPolEq);


% Reachability Analysis ---------------------------------------------------
      
tic
R = reach(vanderPol, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% simulation settings
simOpt.points = 60;
% random simulation
simRes = simulateRandom(vanderPol, params, simOpt);


% Visualization -----------------------------------------------------------

projDim = [1 2];
figure; hold on; box on;

% plot reachable sets
useCORAcolors("CORA:contDynamics")
plot(R,projDim,'DisplayName','Reachable set');

% plot initial set
plot(R(1).R0,projDim, 'DisplayName','Initial set');

% plot simulation results      
plot(simRes,projDim,'DisplayName','Simulations');

% label plot
xlabel(['x_{',num2str(projDim(1)),'}']);
ylabel(['x_{',num2str(projDim(2)),'}']);


% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
