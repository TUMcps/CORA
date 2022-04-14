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
%    no
%
% Outputs:
%    completed - boolean 
%
% References:
%    [1] M. Althoff and J. M. Dolan. Online verification of automated
%        road vehicles using reachability analysis.
%        IEEE Transactions on Robotics, 30(4):903-918, 2014.

% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  23-April-2020 (restucture params/options)
% Last revision:---

%------------- BEGIN CODE --------------

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

vehicle = nonlinearSys(@vmodel_A_bicycle_linear_controlled);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(vehicle, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% simulation settings
simOpt.points = 60;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 400;

% random simulation
simRes = simulateRandom(vehicle, params, simOpt);


% Visualization -----------------------------------------------------------

dims = {[1 2],[3 4],[5 6]};
ref = {[17 18],[19 20],[21 22]};

for k = 1:length(dims)
    
    figure; hold on; box on
    projDims = dims{k}; projRef = ref{k};

    % plot reachable sets 
    plot(R,projDims,'FaceColor',[.8 .8 .8],'EdgeColor','none','Order',3);
    
    % plot initial set
    plot(params.R0,projDims,'w','Filled',true,'EdgeColor','k');
    
    % plot simulation results     
    plot(simRes,projDims);
    
    % plot reference trajectory
    plot(params.u(projRef(1),:),params.u(projRef(2),:),'r','LineWidth',2);

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
completed = 1;

%------------- END OF CODE --------------