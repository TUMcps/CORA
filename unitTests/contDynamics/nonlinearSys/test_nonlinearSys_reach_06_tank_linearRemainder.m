function res = test_nonlinearSys_reach_06_tank_linearRemainder
% test_nonlinearSys_reach_06_tank_linearRemainder - example of
%    nonlinear reachability  analysis;
%
% This example can be found in [1, Sec. 3.4.5] or in [2].
%
% Syntax:
%    test_nonlinearSys_reach_06_tank_linearRemainder
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Althoff. Reachability analysis and its application to the safety 
%        assessment of autonomous cars", Dissertation, TUM 2010.
%    [2] M. Althoff, O. Stursberg, and M. Buss. Reachability analysis of
%        nonlinear systems with uncertain parameters using conservative
%        linearization. In Proc. of the 47th IEEE Conference on
%        Decision and Control, pages 4042–4048, 2008.

% Authors:       Victor Gassmann
% Written:       23-May-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dim_x=6;
params.tFinal=20; %final time
params.R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]);
params.U = zonotope([0,0.005]);


% Reachability Settings ---------------------------------------------------

options.timeStep=4; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=50; %zonotope order
options.reductionInterval=1e3;
options.maxError = ones(dim_x,1);

options.alg = 'lin';
options.tensorOrder = 2;

% System Dynamics----------------------------------------------------------

tank = nonlinearSys(@tank6Eq); %initialize tank system


% Reachability Analysis ---------------------------------------------------

tx1 = tic;
R_wo_linear = reach(tank, params, options); %with normal remainder
tComp1 = toc(tx1);
disp(['computation time of reachable set with normal lagrange remainder: ',num2str(tComp1)]);
tx2 = tic;
options.alg = 'linRem';
options.intermediateOrder=5;
R = reach(tank, params, options); %remainder added to system matrices
tComp2 = toc(tx2);
disp(['computation time of reachable set with remainder added to system matrix: ',num2str(tComp2)]);


% Simulation --------------------------------------------------------------

simRes = simulateRandom(tank, params);


% Visualization -----------------------------------------------------------

plotting = false;

if plotting

    for plotRun=1:3
        % plot different projections
        if plotRun==1
            projDims=[1 2];
        elseif plotRun==2
            projDims=[3 4];    
        elseif plotRun==3
            projDims=[5 6]; 
        end 

        figure; hold on; box on;

        % plot reachable set (normal lagrange remainder)
        plot(R_wo_linear,projDims,'b');

        % plot reachable sets (lagrange remainder added to system matrices
        % (A,B))
        plot(R,projDims,'r');


        %plot initial set
        plot(params.R0,projDims,'k','FaceColor','w');


        %plot simulation results      
        plot(simRes,projDims,'k');

        %label plot
        xlabel(['x_{',num2str(projDims(1)),'}']);
        ylabel(['x_{',num2str(projDims(2)),'}']);
    end

end

% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
