function completed = example_nonlinear_reach_05_autonomousCar()
% example_nonlinear_reach_05_autonomousCar - example of 
%     nonlinear reachability analysis for following a reference trajectory; 
%     this example is also a unit test function.
%
%     This example is similar to the one in [1].
%
%     One of the differences is that computational accelerations such as 
%     parallelization are not activated for this example.
%     Further accelerations, such as taking advantage of monotonicity
%     in the Lagrange remainder are also not considered.
%
% Syntax:
%    completed = example_nonlinear_reach_05_autonomousCar()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff and J. M. Dolan. Online verification of automated
%        road vehicles using reachability analysis.
%        IEEE Transactions on Robotics, 30(4):903-918, 2014.

% Authors:       Matthias Althoff
% Written:       18-August-2016
% Last update:   23-April-2020 (restucture params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 4.01;

% intial set
params.R0 = zonotope([[0; 0; 0; 22; 0 ; 0; -2.1854; 0],...
                      0.05*diag([1, 1, 1, 1, 1, 1, 1, 1])]);
                  
% reference trajectory
params.u = uTRansVec4CASreach();

% uncertain inputs
params.U = zonotope([0*params.u(:,1), 0.05*diag([ones(5,1);zeros(21,1)])]);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.01;
options.taylorTerms = 5;
options.zonotopeOrder = 200;

options.alg = 'lin';
options.tensorOrder = 2;


% System Dynamics ---------------------------------------------------------

vehicle = nonlinearSys(@vmodel_A_bicycle_linear_controlled,8,26);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(vehicle, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% simulation settings
simOpt.points = 10;
% random simulation
simRes = simulateRandom(vehicle, params, simOpt);


% Visualization -----------------------------------------------------------

dims = {[1 2],[3 4],[5 6]};
ref = {[17 18],[19 20],[21 22]};

for k = 1:length(dims)
    
    figure; hold on; box on
    projDim = dims{k}; projRef = ref{k};

    % plot reachable sets
    useCORAcolors("CORA:contDynamics")
    plot(R,projDim,'DisplayName','Reachable set');
    
    % plot initial set
    plot(R(1).R0,projDim, 'DisplayName','Initial set');
    
    % plot simulation results      
    plot(simRes,projDim,'DisplayName','Simulations');
    
    % plot reference trajectory
    plot(params.u(projRef(1),:),params.u(projRef(2),:),'r', ...
        'LineWidth',2,'DisplayName','Reference trajectory');

    % label plot
    xlabel(['x_{',num2str(projDim(1)),'}']);
    ylabel(['x_{',num2str(projDim(2)),'}']);
    legend()
end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
