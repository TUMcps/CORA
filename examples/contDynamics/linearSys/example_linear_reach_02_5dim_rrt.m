function res = example_linear_reach_02_5dim_rrt()
% example_linear_reach_02_5dim_rrt - example of linear reachability 
%    analysis with uncertain inputs and simulation with Rapidly Exploring
%    Random Trees, can be found in [1, Sec. 3.2.3].
%
%    The difference to example_linear_reach_01_5dim is that an RRT is used
%    to generate possible behaviors.
%
% Syntax:
%    res = example_linear_reach_02_5dim_rrt()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010

% Authors:       Matthias Althoff
% Written:       23-September-2016
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

params.tFinal = 1;
params.R0 = zonotope([ones(5,1),0.1*eye(5)]);
params.U = 0.5*zonotope([[1;0;0;0.5;-0.5],diag([0.2,0.5,0.2,0.5,0.5])]);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.04; 
options.taylorTerms = 4;
options.zonotopeOrder = 200;


% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;

fiveDimSys = linearSys('fiveDimSys',A,B);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(fiveDimSys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% simulation options
simOpt.points = 100;
simOpt.vertSamp = true;
simOpt.stretchFac = 1;
simOpt.type = 'rrt';
simOpt.R = R;

% simulate using Rapidly Exploring Random Trees
simRes = simulateRandom(fiveDimSys,params,simOpt);


% Visualization -----------------------------------------------------------

figure;
% plot different projections
dims = {[1,2],[3 4]};

for k = 1:length(dims)
    
    subplot(1,2,k); hold on; box on
    useCORAcolors("CORA:contDynamics");
    projDims = dims{k};

    % plot reachable sets 
    plot(R,projDims, 'DisplayName', 'Reachable set');
    
    % plot initial set
    plot(R.R0,projDims,'DisplayName', 'Initial set');
    
    % plot simulation results
    plot(simRes,projDims, 'DisplayName', 'Simulations');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
