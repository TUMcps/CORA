function res = testLongDuration_nonlinearSys_reach_07_VDP_linearRemainder
% testLongDuration_nonlinearSys_reach_07_vanDerPol_linearRemainder - example of
%    nonlinear reachability analysis;
%
% Syntax:  
%    testLongDuration_nonlinearSys_reach_07_vanDerPol_linearRemainder
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%    -

% Author:       Victor Gassmann
% Written:      22-May-2019
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameters --------------------------------------------------------------

params.tFinal=2.5; %final time
Z0{1}=zonotope([1.4;2.3],diag([0.3,0.05])); %initial state for reachability analysis
params.R0 = Z0{1};
params.U = zonotope(0); %input for reachability analysis


% Reachability Settings ---------------------------------------------------

options.timeStep=0.02; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.intermediateOrder = 10;
options.errorOrder = 5;

options.tensorOrder=2;
options.maxError=0.5*[1; 1];
options.reductionInterval=100;


% System Dynamics ---------------------------------------------------------

vanderPol=nonlinearSys(@vanderPolEq); %initialize van-der-Pol oscillator


% Reachability Analysis ---------------------------------------------------

options.alg = 'lin';
tx1 = tic;
%compute reachable set
R_wo_linear = reach(vanderPol,params,options);
tComp1 = toc(tx1);
disp(['computation time of reachable set with normal remainder: ',num2str(tComp1)]);

tx2 = tic;
%compute reachable set
options.alg = 'linRem';
R = reach(vanderPol,params,options);
tComp2 = toc(tx2);
disp(['computation time of reachable set with remainder added to system matrices: ',num2str(tComp2)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;
simRes = simulateRandom(vanderPol, params, simOpt);


% example completed
res = true;

%------------- END OF CODE --------------