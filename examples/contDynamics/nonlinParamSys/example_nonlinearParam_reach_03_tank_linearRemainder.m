function completed = example_nonlinearParam_reach_03_tank_linearRemainder()
% example_nonlinearParam_reach_03_tank_linearRemainder - example of
%     nonlinear reachability analysis with uncertain parameters, taken from
%     [1, Sec. 3.4.5] or in [2].
%
% Syntax:
%    completed = example_nonlinearParam_reach_03_tank_linearRemainder
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, "Reachability analysis and its application to the safety 
%        assessment of autonomous cars", Dissertation, TUM, 2010, 
%    [2] M. Althoff, O. Stursberg, and M. Buss. Reachability analysis
%        of nonlinear systems with uncertain parameters using
%        conservative linearization. CDC 2008

% Authors:       Victor Gassmann
% Written:       23-May-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

dim_x=6;
params.R0=zonotope([[2; 4; 4; 2; 10; 4],1*eye(dim_x)]);
params.U=zonotope([0,0.005]); %input for reachability analysis
params.tFinal=20; %final time
params.paramInt=interval(0.015,0.015); %parameter intervals for reachability analysis


% Reachability Settings ---------------------------------------------------

options.timeStep=1;
options.taylorTerms=4; %number of taylor terms for reachable sets
options.intermediateTerms = 4;
options.zonotopeOrder=10; %zonotope order
options.maxError = 0.1*ones(dim_x,1);
options.reductionInterval=1e3;
options.tensorOrder = 2;
options.alg = 'lin';


% System Dynamics ---------------------------------------------------------

tank = nonlinearSys(@tank6Eq);
tankParam = nonlinParamSys(@tank6paramEq); %with uncertain parameters


% Reachability Analysis ---------------------------------------------------

tx1 = tic;
R_wo_linear = reach(tankParam, params, options); %with normal remainder
tComp1 = toc(tx1);
disp(['computation time of reachable set with normal lagrange remainder: ',num2str(tComp1)]);

tx2 = tic;
options.alg='linRem';
R_param = reach(tankParam, params, options); %remainder added to system matrices
tComp2 = toc(tx2);
disp(['computation time of reachable set with remainder added to system matrix: ',num2str(tComp2)]);


% Simulation --------------------------------------------------------------

simOpt.points = 60;
params = rmfield(params,'paramInt');

simRes = simulateRandom(tank, params, simOpt);


% Visualization -----------------------------------------------------------

plotOrder = 8;
for plotRun=1:3
    % plot different projections
    if plotRun==1
        projectedDims=[1 2];
    elseif plotRun==2
        projectedDims=[3 4];    
    elseif plotRun==3
        projectedDims=[5 6]; 
    end 

    figure; hold on; box on;
    useCORAcolors("CORA:contDynamics", 2)
    
    % reachable set: normal lagrange remainder
    plot(R_wo_linear,projectedDims,'Order',plotOrder);
    
    % reachable set: lagrange remainder added to system matrices (A,B)
    plot(R_param,projectedDims,'Order',plotOrder);
    
    %plot initial set
    plot(R_param.R0,projectedDims);
  
    %plot simulation results      
    plot(simRes,projectedDims);

    %label plot
    xlabel(['x_{',num2str(projectedDims(1)),'}']);
    ylabel(['x_{',num2str(projectedDims(2)),'}']);
end


%example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
