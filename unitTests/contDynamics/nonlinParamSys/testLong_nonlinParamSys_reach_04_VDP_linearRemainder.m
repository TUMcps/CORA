function completed = testLong_nonlinParamSys_reach_04_VDP_linearRemainder
% testLong_nonlinParamSys_reach_04_VDP_linearRemainder - example of
%    nonlinear reachability analysis with uncertain parameters; 
%
% Syntax:
%    completed = testLong_nonlinParamSys_reach_04_VDP_linearRemainder
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Victor Gassmann
% Written:       17-May-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

dim_x=2;
params.tFinal=2.5; %final time
params.R0=zonotope([[1;1],0.1*[0.3 0;0 0.05]]);
params.U=zonotope([0,0.005]);

params.paramInt=interval([0.995;1],[1;1.005]); %parameter intervals for reachability analysis


% Reachability Analysis ---------------------------------------------------

options.timeStep=0.1; %time step size for reachable set computation
options.taylorTerms=5; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.intermediateTerms = 4;
options.maxError = 2*ones(dim_x,1);
options.reductionInterval=1e3;
options.tensorOrder = 2;
options.alg = 'lin';

% System Dynamics ---------------------------------------------------------

vanderPol = nonlinearSys(@vanderPolEq);
vanderPolParam = nonlinParamSys(@vanderPolparamEq); %with uncertain parameters


% Reachability Analysis ---------------------------------------------------

tx1 = tic;
R_wo_linear = reach(vanderPolParam,params, options); %with normal remainder
tComp1 = toc(tx1);
disp(['computation time of reachable set with normal lagrange remainder: ',num2str(tComp1)]);

tx2 = tic;
options.alg='linRem';
R_Param = reach(vanderPolParam,params,options); %remainder added to system matrices
tComp2 = toc(tx2);
disp(['computation time of reachable set with remainder added to system matrix: ',num2str(tComp2)]);


% Simulation --------------------------------------------------------------

params = rmfield(params,'paramInt');
simOpt.points = 60;
simRes = simulateRandom(vanderPol, params, simOpt);


% Visualization -----------------------------------------------------------

% don't plot in suite
plotting = false;

if plotting

    projectedDims = [1 2];
    plotOrder = 20;

    figure; hold on; box on;

    % reachable set: standard lagrange remainder
    plot(R_wo_linear,projectedDims,'b','Order',plotOrder);

    % reachable set: lagrange remainder added to system matrices (A,B)
    plot(R_Param,projectedDims,'r','Order',plotOrder);

    %plot initial set
    plot(params.R0,projectedDims,'k','FaceColor','w');

    %plot simulation results
    plot(simRes,projectedDims,'k');

    %label plot
    xlabel(['x_{',num2str(projectedDims(1)),'}']);
    ylabel(['x_{',num2str(projectedDims(2)),'}']);

end

%example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
