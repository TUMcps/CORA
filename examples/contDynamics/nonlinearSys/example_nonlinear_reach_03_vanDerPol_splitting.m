function completed = example_nonlinear_reach_03_vanDerPol_splitting()
% example_nonlinear_reach_03_vanDerPol - example of nonlinear reachability 
%                                        analysis with splitting
%
% This example can be found in  [1, Sec. 3.4.5] or in [2].
% A new technique for computing this example with less spliiting has been
% published in [3].
%
% Syntax:  
%    completed = example_nonlinear_reach_03_vanDerPol()
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean 
%
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010
%    [2] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008
%    [3] M. Althoff. "Reachability analysis of nonlinear systems using 
%        conservative polynomialization and non-convex sets", HSCC 2013

% Author:       Matthias Althoff
% Written:      26-June-2009
% Last update:  23-April-2020 (restucture params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameters --------------------------------------------------------------

params.tFinal = 3.2;

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
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;

% random simulation
simRes = simulateRandom(vanderPol, params, simOpt);


% Visualization -----------------------------------------------------------

dim = [1 2];
plotOrder = 20;
    
figure; hold on;

% plot reachable sets 
plot(R,dim,'FaceColor',[.8 .8 .8],'EdgeColor','none','Order',plotOrder);

% plot initial set
plot(params.R0,dim,'w','Filled',true,'EdgeColor','k');

% plot simulation results   
plot(simRes,dim);

% label plot
xlabel(['x_{',num2str(dim(1)),'}']);
ylabel(['x_{',num2str(dim(2)),'}']);


% example completed
completed = 1;

%------------- END OF CODE --------------