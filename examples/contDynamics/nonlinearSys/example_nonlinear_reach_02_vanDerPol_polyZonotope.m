function completed = example_nonlinear_reach_02_vanDerPol_polyZonotope()
% example_nonlinear_reach_02_vanDerPol_polyZonotope - example of non-linear 
%                                                     reachability analysis 
%
%    Example from [1] comparing reachability analysis for non-linear
%    systems using zonotopes and polynomial zonotopes for the Van-der-Pol
%    oscillator system.
%
% Syntax:
%    completed = example_nonlinear_reach_02_vanDerPol_polyZonotope()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References: 
%   [1] N. Kochdumper et al. "Sparse Polynomial Zonotopes: A Novel Set 
%       Representation for Reachability Analysis"

% Authors:       Niklas Kochdumper
% Written:       02-January-2020
% Last update:   23-April-2020 (restucture params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 6.74;                              % final time
R0 = zonotope([[1.4;2.4], diag([0.17, 0.06])]);    % initial set


% Reachability Settings ---------------------------------------------------

% settings
options.timeStep = 0.005;                           
options.taylorTerms = 4;                            
options.zonotopeOrder = 50;       
options.intermediateOrder = 50;
options.errorOrder = 20;

% reachability algorithm
options.alg = 'poly';
options.tensorOrder = 3;

% special settings for polynomial zonotopes
polyZono.maxDepGenOrder = 50;
polyZono.maxPolyZonoRatio = 0.01;
polyZono.restructureTechnique = 'reducePca';


% System Dynamics ---------------------------------------------------------

vanderPol = nonlinearSys(@vanderPolEq);


% Reachability Analysis (zonotope) ----------------------------------------
      
% adapted options
params.R0 = R0;

% compute reachable set
tic
R = reach(vanderPol, params, options);
tComp = toc;
disp(['computation time (zonotope): ',num2str(tComp)]);


% Reachability Analysis (polynomial zonotope) -----------------------------
      
% adapted options
params.R0 = polyZonotope(R0);

options.polyZono = polyZono;

% compute reachable set
tic
Rpoly = reach(vanderPol, params, options);
tComp = toc;
disp(['computation time (polynomial zonotope): ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% simulation settings
simOpt.points = 20;
% random simulation
simRes = simulateRandom(vanderPol, params, simOpt);

% Visualization -----------------------------------------------------------
    
figure; hold on; box on;
projDim = [1 2];

% plot reachable set (zonotope)

useCORAcolors("CORA:contDynamics", 2)
plot(R,projDim,'DisplayName','zonotope');
plot(Rpoly,projDim,'DisplayName','polynomial zonotope');

% plot initial set
plot(R(1).R0,projDim, 'DisplayName','Initial set');

% plot simulation results      
plot(simRes,projDim,'DisplayName','Simulations');

% label plot
xlabel('x_1');
ylabel('x_2');

legend();


% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
