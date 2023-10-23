function res = testLong_nonlinearSys_reach_07_VDP_linearRemainder
% testLong_nonlinearSys_reach_07_VDP_linearRemainder -
%    example of nonlinear reachability analysis;
%
% Syntax:
%    testLong_nonlinearSys_reach_07_VDP_linearRemainder
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Victor Gassmann
% Written:       22-May-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 2.5;
Z0{1} = zonotope([1.4;2.3],diag([0.3,0.05]));
params.R0 = Z0{1};
params.U = zonotope(0);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.02;
options.taylorTerms = 4;
options.zonotopeOrder = 10;

options.tensorOrder = 2;
options.maxError = 0.5*[1; 1];
options.reductionInterval = 100;


% System Dynamics ---------------------------------------------------------

vanderPol = nonlinearSys(@vanderPolEq);


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
options.intermediateOrder = 10;
R = reach(vanderPol,params,options);
tComp2 = toc(tx2);
disp(['computation time of reachable set with remainder added to system matrices: ',num2str(tComp2)]);


% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
