function completed = example_nonlinearParam_reach_04_vanderPol_linearRemainder()
% example_nonlinearParam_reach_04_vanderPol_linearRemainder - example
%    of nonlinear reachability analysis with uncertain parameters; 
%
% Syntax:  
%    example_nonlinearParam_reach_04_vanderPol_linearRemainder
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%

% Author:       Victor Gassmann
% Written:      17-May-2019
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameters --------------------------------------------------------------

dim_x=2;
params.R0=zonotope([[1.4;2.3],0.1*[0.3 0;0 0.05]]);
params.U=zonotope([0,0.005]);
params.tFinal=2.5; %final time


% Reachability Settings ---------------------------------------------------

options.timeStep=0.1; %time step size for reachable set computation
options.taylorTerms=5; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.intermediateOrder = 4;
options.maxError = 1*ones(dim_x,1);
options.reductionInterval=1e3;
options.tensorOrder = 2;
options.verbose = true;
options.alg = 'lin';

%parameter intervals for reachability analysis
options.paramInt=interval([0.995;1],[1;1.005]);


% System Dynamics ---------------------------------------------------------

vanderPol = nonlinearSys(@vanderPolEq);
vanderPolParam = nonlinParamSys(@vanderPolparamEq); %with uncertain parameters


% Reachability Analysis ---------------------------------------------------

tx1 = tic;
R_wo_linear = reach(vanderPolParam, params, options); %with normal remainder
tComp1 = toc(tx1);
disp(['computation time of reachable set with normal lagrange remainder: ',num2str(tComp1)]);

tx2 = tic;
options.alg='linRem';
R_param = reach(vanderPolParam, params, options); %remainder added to system matrices
tComp2 = toc(tx2);
disp(['computation time of reachable set with remainder added to system matrix: ',num2str(tComp2)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;

simRes = simulateRandom(vanderPol, params, simOpt);


% Simulation --------------------------------------------------------------

projDims=[1 2];
plotOrder = 20;

figure; hold on; box on;

% reachable set: normal lagrange remainder
plot(R_wo_linear,projDims,'FaceColor',[.5 .5 .5],'Order',plotOrder,'EdgeColor','none');

% reachable set: lagrange remainder added to system matrices (A,B)
plot(R_param,projDims,'FaceColor',[.8 .8 .8],'Order',plotOrder,'EdgeColor','none');

% plot initial set
plot(params.R0,projDims,'w','Filled',true,'EdgeColor','k');

% plot simulation results
plot(simRes,projDims,'k');

% label plot
xlabel(['x_{',num2str(projDims(1)),'}']);
ylabel(['x_{',num2str(projDims(2)),'}']);


% example completed
completed = 1;

%------------- END OF CODE --------------
