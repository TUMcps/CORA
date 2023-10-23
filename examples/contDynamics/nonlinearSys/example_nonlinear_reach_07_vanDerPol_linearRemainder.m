function completed = example_nonlinear_reach_07_vanDerPol_linearRemainder()
% example_nonlinear_reach_07_vanDerPol_linearRemainder - example of
%    nonlinear reachability analysis;
%
% Syntax:
%    completed = example_nonlinear_reach_07_vanDerPol_linearRemainder()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 

% Authors:       Victor Gassmann
% Written:       22-May-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal=2.5; %final time
Z0{1}=zonotope([1.4;2.3],[0.3 0;0 0.05]);
params.R0 = zonoBundle(Z0);
params.U = zonotope(0);


% Reachability Settings ---------------------------------------------------

options.timeStep=0.02; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order

options.tensorOrder=2;
options.maxError=0.1*[1; 1];
options.reductionInterval=100;
options.verbose = true;


% System Dynamics ---------------------------------------------------------

vanderPol=nonlinearSys(@vanderPolEq); %initialize van-der-Pol oscillator


% Reachability Analysis ---------------------------------------------------
      
tx1 = tic;
%compute reachable set
options.alg = 'lin';
R_wo_linear = reach(vanderPol, params, options);
tComp1 = toc(tx1);
disp(['computation time of reachable set with ' ...
    'normal remainder: ',num2str(tComp1)]);

tx2 = tic;
%compute reachable set
options.alg = 'linRem';
options.intermediateOrder = 10;
R_linRem = reach(vanderPol, params, options);
tComp2 = toc(tx2);
disp(['computation time of reachable set with '...
    'remainder added to system matrices: ',num2str(tComp2)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simRes = simulateRandom(vanderPol, params, simOpt);


% Visualization -----------------------------------------------------------

projDims = [1 2];
    
figure; hold on; box on;

% plot reachable set (normal lagrange remainder)
useCORAcolors("CORA:contDynamics", 2)
plot(R_wo_linear,projDims,'DisplayName','Normal Remainder');

% plot reachable sets (lagrange remainder added to system matrices (A,B))
plot(R_linRem,projDims,'DisplayName','Linear Remainder');

% plot initial set
plot(R_linRem(1).R0,projDims,'DisplayName','Initial set');

% plot simulation results      
plot(simRes,projDims,'DisplayName','Simulations');

%label plot
xlabel(['x_{',num2str(projDims(1)),'}']);
ylabel(['x_{',num2str(projDims(2)),'}']);
legend()

%example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
