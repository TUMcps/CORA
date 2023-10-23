function completed = example_nonlinearDA_reach_04_powerSystem_14bus()
% example_nonlinearDA_reach_04_powerSystem_14bus - example of 
%    nonlinear-differential-algebraic reachability analysis: 14-bus power 
%    system from [1].
%
% Syntax:
%    completed = example_nonlinearDA_reach_04_powerSystem_14bus()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, "Formal and Compositional Analysis of Power Systems 
%        using Reachable Sets", IEEE Transactions on Power Systems 29 (5), 
%        2014, 2270-2280

% Authors:       Matthias Althoff
% Written:       26-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Bus Parameters ----------------------------------------------------------

nrOfGenerators = 5;
nrOfBuses = 14;
nrOfInputs = nrOfGenerators + nrOfBuses;


% Reachability Settings ---------------------------------------------------

options.taylorTerms = 6;
options.zonotopeOrder = 400;
options.tensorOrder = 2;
options.errorOrder = 3;
options.intermediateOrder = 2;
options.lagrangeRem.simplify = 'optimize';

options.maxError_x = [zeros(5,1); 0.6*ones(5,1); zeros(5,1)];
options.maxError_y = 0.01*ones(27,1);


% System Dynamics ---------------------------------------------------------

% set path
path = [CORAROOT filesep 'models' filesep 'powerSystemsConverted'];

% create models (normal and faulty operation)
if ~isfile([path filesep 'IEEE14_model.mat'])
    powerSystem2cora('IEEE14');
end
if ~isfile([path filesep 'IEEE14_fault_model.mat'])
    powerSystem2cora('IEEE14_fault');
end

% load models (normal and faulty operation)
load([path filesep 'IEEE14_model.mat'], 'IEEE14_model');
load([path filesep 'IEEE14_fault_model.mat'], 'IEEE14_fault_model');


% Reachability Analysis ---------------------------------------------------

%% normal operation
% time horizon and time step
params.tFinal = 0.1; 
options.timeStep = 0.005;
timeSteps = params.tFinal/options.timeStep;
% input
u0_load = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
u0_gen = [2; 0.4; 0; 0; 0];
params.u = [u0_gen; u0_load]*ones(1,timeSteps);
% parameters
x0guess = [0*ones(nrOfGenerators,1); 377*ones(nrOfGenerators,1); 1*ones(nrOfGenerators,1)];
y0guess = [ones(nrOfBuses, 1); zeros(nrOfBuses-1, 1)];

[x0,params.y0guess] = steadyState(IEEE14_model, x0guess, y0guess, params.u);

R0 = zonotope([x0, diag([1e-2*ones(nrOfGenerators,1); ...
    1e-1*ones(nrOfGenerators,1); 1e-3*ones(nrOfGenerators,1)])]);
params.R0 = R0;

params.U = zonotope(zeros(nrOfInputs,1));


% computation of reachable set
tic
R_1 = reach(IEEE14_model, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%% operation under fault
% time horizon and time step
params.tStart = 0.1; 
params.tFinal = 0.13; 
options.timeStep = 0.001; 
timeSteps = (params.tFinal - params.tStart)/options.timeStep;
% input
u0_gen = [0; 0.4; 0; 0; 0];
params.u = [u0_gen; u0_load]*ones(1,timeSteps);
% initial set is final set of normal operation
params.R0 = R_1.timePoint.set{end}; %initial state for reachability analysis
% computation of reachable set
tic
R_2 = reach(IEEE14_fault_model, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%% operation after fixing the fault
% time horizon and time step
params.tStart = 0.13; 
params.tFinal = 2; 
options.timeStep = 0.005; 
timeSteps = (params.tFinal - params.tStart)/options.timeStep;
% input
u0_gen = [2; 0.4; 0; 0; 0];
params.u = [u0_gen; u0_load]*ones(1,timeSteps);
% initial set is final set of normal operation
params.R0 = R_2.timePoint.set{end}; %initial state for reachability analysis
% computation of reachable set
tic
R_3 = reach(IEEE14_model, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

nrRuns = 10;
for iRun = 1:nrRuns
    simOpt.points = 1;
    
    % normal operation
    params.R0 = R0;
    params.tStart = 0;
    params.tFinal = 0.1; 
    u0_load = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
    u0_gen = [2; 0.4; 0; 0; 0];
    params.u = [u0_gen; u0_load];
    simRes_1{iRun} = simulateRandom(IEEE14_model, params, simOpt);

    % operation under fault
    params.tStart = 0.1; 
    params.tFinal = 0.13; 
    u0_gen = [0; 0.4; 0; 0; 0];
    params.u = [u0_gen; u0_load];
    % initial set is final solution of normal operation
    params.R0 = zonotope(simRes_1{iRun}.x{1}(end,:)');
    simRes_2{iRun} = simulateRandom(IEEE14_fault_model, params, simOpt);
    
    % operation after fixing the fault
    params.tStart = 0.13; 
    params.tFinal = 2; 
    u0_gen = [2; 0.4; 0; 0; 0];
    params.u = [u0_gen; u0_load];
    % initial set is final solution of normal operation
    params.R0 = zonotope(simRes_2{iRun}.x{1}(end,:)');
    simRes_3{iRun} = simulateRandom(IEEE14_model, params, simOpt);
end
    
% Visualization -----------------------------------------------------------

doPlot = false;
if ~doPlot
    disp('Visualiation is turned off for this example for time reasons. You can enable it by setting the ''doPlot'' variable to true.')
    return
end

projDims = {[1 6],[2 7],[3 8],[4 9],[5 10],[11 12],[13 14],[6 11],...
    [7 12],[8 13],[9 14],[1 11],[2 12],[3 13],[4 14],[1 2],[3 4],[5 6],[7 8]};

%plot dynamic variables
for plotRun=1:19
    projDim = projDims{plotRun};

    figure; hold on; box on;

    % plot reachable sets
    plot(R_1, projDim);
    plot(R_2, projDim);
    plot(R_3, projDim);
    
    % plot initial set
    plot(params.R0,projDim,'EdgeColor','k','FaceColor','w');

    % plot simulation results
    for j=1:nrRuns
        plot(simRes_1{j},projDim);
        plot(simRes_2{j},projDim);
        plot(simRes_3{j},projDim);
    end

    % label plot
    xlabel(['x_{',num2str(projDim(1)),'}']);
    ylabel(['x_{',num2str(projDim(2)),'}']);
end

projDims = {[1 15],[2 16],[3 17],[4 18],[5 19],[6 20],[7 21],[8 22],...
    [9 23],[10 24],[11 25],[12 26],[13 27],[14 27],[1 2],[3 4],[5 6],[7 8]};

%plot constraint variables
for plotRun=1:18
    projDim = projDims{plotRun};
    
    figure; hold on; box on;

    % plot reachable sets
    plot(R_1,projDim,'Set','y');
    plot(R_2,projDim,'Set','y');
    plot(R_3,projDim,'Set','y');

    % plot initial set
%     plot(params.R0,projDim,'EdgeColor','k','FaceColor','w');

    % plot simulation results
    for j=1:nrRuns
        plot(simRes_1{j},projDim,'Traj','a');
        plot(simRes_2{j},projDim,'Traj','a');
        plot(simRes_3{j},projDim,'Traj','a');
    end

    % label plot
    xlabel(['y_{',num2str(projDim(1)),'}']);
    ylabel(['y_{',num2str(projDim(2)),'}']);

end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
