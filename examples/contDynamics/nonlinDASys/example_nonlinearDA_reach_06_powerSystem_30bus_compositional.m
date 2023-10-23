function completed = example_nonlinearDA_reach_06_powerSystem_30bus_compositional()
% example_nonlinearDA_reach_06_powerSystem_30bus_compositional - example of 
%    nonlinear-differential-algebraic reachability analysis: 30-bus power 
%    system from [1], where the abstraction error is computed 
%    compositionally.
%
% Syntax:
%    completed = example_nonlinearDA_reach_06_powerSystem_30bus_compositional()
%
% Inputs:
%    ---
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, "Formal and Compositional Analysis of Power Systems 
%        using Reachable Sets", IEEE Transactions on Power Systems 29 (5), 
%        2014, 2270-2280

% Authors:       Matthias Althoff
% Written:       09-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Bus Parameters ----------------------------------------------------------

nrOfGenerators = 6;
nrOfBuses = 30;
nrOfInputs = nrOfGenerators + nrOfBuses;


% Reachability Settings ---------------------------------------------------

options.taylorTerms = 6;
options.zonotopeOrder = 400;
options.tensorOrder = 2;
options.errorOrder = 5;
options.intermediateOrder = 2;
options.lagrangeRem.simplify = 'optimize';

options.maxError_x = [zeros(nrOfGenerators,1); 0.4*ones(nrOfGenerators,1); zeros(nrOfGenerators,1)]; % for comparison reasons
options.maxError_y = 0.04*ones(2*nrOfBuses-1,1);


% System Dynamics ---------------------------------------------------------

% set path
path = [CORAROOT filesep 'models' filesep 'powerSystemsConverted'];

% create models (normal and faulty operation)
if ~isfile([path filesep 'IEEE30_model.mat'])
    powerSystem2cora('IEEE30');
    load([path filesep 'IEEE30_model'], 'IEEE30_model');
    derivatives(IEEE30_model, options); % derivatives are not automatically created
end
if ~isfile([path filesep 'IEEE30_fourParts_sub1_model.mat'])
    powerSystem2cora('IEEE30_fourParts_sub1');
    load([path filesep 'IEEE30_fourParts_sub1_model'], 'IEEE30_fourParts_sub1_model');
    derivatives(IEEE30_fourParts_sub1_model, options); % derivatives are not automatically created for subsystems
end
if ~isfile([path filesep 'IEEE30_fourParts_sub2_model.mat'])
    powerSystem2cora('IEEE30_fourParts_sub2');
    load([path filesep 'IEEE30_fourParts_sub2_model'], 'IEEE30_fourParts_sub2_model');
    derivatives(IEEE30_fourParts_sub2_model, options); % derivatives are not automatically created for subsystems
end
if ~isfile([path filesep 'IEEE30_fourParts_sub3_model.mat'])
    powerSystem2cora('IEEE30_fourParts_sub3');
    load([path filesep 'IEEE30_fourParts_sub3_model'], 'IEEE30_fourParts_sub3_model');
    derivatives(IEEE30_fourParts_sub3_model, options); % derivatives are not automatically created for subsystems
end
if ~isfile([path filesep 'IEEE30_fourParts_sub4_model.mat'])
    powerSystem2cora('IEEE30_fourParts_sub4');
    load([path filesep 'IEEE30_fourParts_sub4_model'], 'IEEE30_fourParts_sub4_model');
    derivatives(IEEE30_fourParts_sub4_model, options); % derivatives are not automatically created for subsystems
end
if ~isfile([path filesep 'IEEE30_fault_model.mat'])
    powerSystem2cora('IEEE30_fault');
    load([path filesep 'IEEE30_fault_model'], 'IEEE30_fault_model');
    derivatives(IEEE30_fault_model, options); % derivatives are not automatically created
end
if ~isfile([path filesep 'IEEE30_fault_fourParts_sub1_model.mat'])
    powerSystem2cora('IEEE30_fault_fourParts_sub1');
    load([path filesep 'IEEE30_fault_fourParts_sub1_model'], 'IEEE30_fault_fourParts_sub1_model');
    derivatives(IEEE30_fault_fourParts_sub1_model, options); % derivatives are not automatically created for subsystems
end

% load models (normal and faulty operation)
load([path filesep 'IEEE30_model'], 'IEEE30_model');
load([path filesep 'IEEE30_fourParts_sub1_model'], 'IEEE30_fourParts_sub1_model');
load([path filesep 'IEEE30_fourParts_sub2_model'], 'IEEE30_fourParts_sub2_model');
load([path filesep 'IEEE30_fourParts_sub3_model'], 'IEEE30_fourParts_sub3_model');
load([path filesep 'IEEE30_fourParts_sub4_model'], 'IEEE30_fourParts_sub4_model');
load([path filesep 'IEEE30_fault_model'], 'IEEE30_fault_model');
load([path filesep 'IEEE30_fault_fourParts_sub1_model'], 'IEEE30_fault_fourParts_sub1_model');

% rewrite systems as cells
sys{1} = IEEE30_model;
sys{2} = IEEE30_fourParts_sub1_model;
sys{3} = IEEE30_fourParts_sub2_model;
sys{4} = IEEE30_fourParts_sub3_model;
sys{5} = IEEE30_fourParts_sub4_model;
sys_fault{1} = IEEE30_fault_model;
sys_fault{2} = IEEE30_fault_fourParts_sub1_model;
sys_fault{3} = sys{3}; % this subsystem has no fault
sys_fault{4} = sys{4}; % this subsystem has no fault
sys_fault{5} = sys{5}; % this subsystem has no fault

% build indices for the input, state, and constraints relating the full 
% full system
load IEEE30.mat bus
busFull = bus;
% subsystem 1
load IEEE30_fourParts_sub1.mat bus
bus_subsystem{1} = bus;
% subsystem 2
load IEEE30_fourParts_sub2.mat bus
bus_subsystem{2} = bus;
% subsystem 3
load IEEE30_fourParts_sub3.mat bus
bus_subsystem{3} = bus;
% subsystem 4
load IEEE30_fourParts_sub4.mat bus
bus_subsystem{4} = bus;
% create cell of subsystems
subsystem(1:4) = sys(2:5);
% compute indices for subsystems
index = indexForSubsystems(sys{1}, subsystem, busFull, bus_subsystem, 3);

% faulty system: build indices for the input, state, and constraints 
% relating the full 
% full system
load IEEE30_fault.mat bus
busFull = bus;
% subsystem 1
load IEEE30_fault_fourParts_sub1.mat bus
bus_subsystem{1} = bus;
% subsystem 2; this subsystem has no fault
load IEEE30_fourParts_sub2.mat bus
bus_subsystem{2} = bus;
% subsystem 3; this subsystem has no fault
load IEEE30_fourParts_sub3.mat bus
bus_subsystem{3} = bus;
% subsystem 4; this subsystem has no fault
load IEEE30_fourParts_sub4.mat bus
bus_subsystem{4} = bus;
% create cell of subsystems
subsystem(1:4) = sys_fault(2:5);
% compute indices for subsystems
index_fault = indexForSubsystems(sys_fault{1}, subsystem, busFull, bus_subsystem, 3);

% Reachability Analysis ---------------------------------------------------

%% normal operation
% time horizon and time step
params.tFinal = 0.1; 
options.timeStep = 0.005;
timeSteps = params.tFinal/options.timeStep;
% provide index for compositional computation
options.index = index;
% load model parameters to obtain Vgen
load IEEE30.mat VM bus
options.Vgen = VM(bus.generator);
% store Jacobian and Hessians of subsystems
for iSys = 2:5
    sys{iSys}.setHessian('int');
    options.subsystemHessian{iSys-1} = sys{iSys}.hessian;
    options.subsystemJacobian{iSys-1} = sys{iSys}.jacobian;
end
% input
u0_load = zeros(nrOfBuses,1);
u0_gen = [2.6; 0.4; 0; 0; 0; 0];
params.u = [u0_gen; u0_load]*ones(1,timeSteps);
% parameters
x0guess = [0*ones(nrOfGenerators,1); 377*ones(nrOfGenerators,1); 1*ones(nrOfGenerators,1)];
y0guess = [ones(nrOfBuses, 1); zeros(nrOfBuses-1, 1)];
 
[x0,params.y0guess] = steadyState(IEEE30_model, x0guess, y0guess, params.u);
 
R0 = zonotope([x0, diag([5e-3*ones(nrOfGenerators,1); ...
     1e-1*ones(nrOfGenerators,1); 1e-3*ones(nrOfGenerators,1)])]);
params.R0 = R0;

params.U = zonotope(zeros(nrOfInputs,1));


% computation of reachable set
tic
R_1 = reach(IEEE30_model, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%% operation under fault
% time horizon and time step
params.tStart = 0.1; 
params.tFinal = 0.13; 
%params.tFinal = 0.15; 
options.timeStep = 0.001; 
timeSteps = (params.tFinal - params.tStart)/options.timeStep;
% provide index for compositional computation
options.index = index_fault;
% load model parameters to obtain Vgen
load IEEE30_fault.mat VM bus
options.Vgen = VM(bus.generator);
% store Jacobian and Hessians of subsystems
sys_fault{2}.setHessian('int');
sys_fault{3}.setHessian('int');
options.subsystemHessian{1} = sys_fault{2}.hessian;
options.subsystemHessian{2} = sys_fault{3}.hessian;
options.subsystemJacobian{1} = sys_fault{2}.jacobian;
options.subsystemJacobian{2} = sys_fault{3}.jacobian;
% input
u0_gen = [0; 0.4; 0; 0; 0; 0];
params.u = [u0_gen; u0_load]*ones(1,timeSteps);
% initial set is final set of normal operation
params.R0 = R_1.timePoint.set{end}; %initial state for reachability analysis
% computation of reachable set
tic
R_2 = reach(IEEE30_fault_model, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%% operation after fixing the fault
% time horizon and time step
params.tStart = 0.13;
%params.tStart = 0.15;
params.tFinal = 2; 
options.timeStep = 0.005; 
timeSteps = round((params.tFinal - params.tStart)/options.timeStep);
% provide index for compositional computation
options.index = index;
% load model parameters to obtain Vgen
load IEEE30.mat VM bus
options.Vgen = VM(bus.generator);
% store Jacobian and Hessians of subsystems
options.subsystemHessian{1} = sys{2}.hessian;
options.subsystemHessian{2} = sys{3}.hessian;
options.subsystemJacobian{1} = sys{2}.jacobian;
options.subsystemJacobian{2} = sys{3}.jacobian;
% input
u0_gen = [2.6; 0.4; 0; 0; 0; 0];
params.u = [u0_gen; u0_load]*ones(1,timeSteps);
% initial set is final set of normal operation
params.R0 = R_2.timePoint.set{end}; %initial state for reachability analysis
% computation of reachable set
tic
R_3 = reach(IEEE30_model, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%% operation after fixing the fault (larger time step size)
% time horizon and time step
params.tStart = 2; 
params.tFinal = 5; 
options.timeStep = 0.02; 
timeSteps = round((params.tFinal - params.tStart)/options.timeStep);
% input
params.u = [u0_gen; u0_load]*ones(1,timeSteps);
% initial set is final set of normal operation
params.R0 = R_3.timePoint.set{end}; %initial state for reachability analysis
% computation of reachable set
tic
R_4 = reach(IEEE30_model, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

save IEEE30reach

% Simulation --------------------------------------------------------------

nrRuns = 10;
for iRun = 1:nrRuns
    simOpt.points = 1;
    
    % normal operation
    params.R0 = R0;
    params.tStart = 0;
    params.tFinal = 0.1; 
    u0_gen = [2.6; 0.4; 0; 0; 0; 0];
    params.u = [u0_gen; u0_load];
    simRes_1{iRun} = simulateRandom(IEEE30_model, params, simOpt);

    % operation under fault
    params.tStart = 0.1; 
    params.tFinal = 0.13; 
    %params.tFinal = 0.15;
    u0_gen = [0; 0.4; 0; 0; 0; 0];
    params.u = [u0_gen; u0_load];
    % initial set is final solution of normal operation
    params.R0 = zonotope(simRes_1{iRun}.x{1}(end,:)');
    simRes_2{iRun} = simulateRandom(IEEE30_fault_model, params, simOpt);
    
    % operation after fixing the fault
    params.tStart = 0.13; 
    %params.tFinal = 0.15;
    params.tFinal = 5; 
    u0_gen = [2.6; 0.4; 0; 0; 0; 0];
    params.u = [u0_gen; u0_load];
    % initial set is final solution of normal operation
    params.R0 = zonotope(simRes_2{iRun}.x{1}(end,:)');
    simRes_3{iRun} = simulateRandom(IEEE30_model, params, simOpt);
end
save IEEE30sim
% Visualization -----------------------------------------------------------

doPlot = false;
if ~doPlot
    disp('Visualiation is turned off for this example for time reasons. You can enable it by setting the ''doPlot'' variable to true.')
    return
end

%plot dynamic variables
projDimDyn{1} = [1 7];
projDimDyn{2} = [2 8];
projDimDyn{3} = [3 9];
projDimDyn{4} = [4 10];
projDimDyn{5} = [5 11];
projDimDyn{6} = [12 13];
projDimDyn{7} = [14 15];      
projDimDyn{8} = [16 17];
projDimDyn{9} = [11 18];
projDimDyn{10} = [8 14];
projDimDyn{11} = [9 15];      
projDimDyn{12} = [1 13];
projDimDyn{13} = [2 14];
projDimDyn{14} = [3 15];      
projDimDyn{15} = [4 16];   
projDimDyn{16} = [1 2];
projDimDyn{17} = [3 4];
projDimDyn{18} = [5 6];
projDimDyn{19} = [7 8];      


for plotRun=1:length(projDimDyn)

    h_dyn{plotRun} = figure; hold on; box on;

    % plot reachable sets
    plot(R_4, projDimDyn{plotRun}, 'Order', 10, 'FaceColor',[.6 .6 .6]);
    plot(R_3, projDimDyn{plotRun}, 'Order', 10, 'FaceColor',[.6 .6 .6]);
    plot(R_2, projDimDyn{plotRun}, 'Order', 10, 'FaceColor',[.8 .8 .8]);
    plot(R_1, projDimDyn{plotRun}, 'Order', 10, 'FaceColor',[.6 .6 .6]);
    
    % plot initial set
    plot(R0,projDimDyn{plotRun},'EdgeColor','k','FaceColor','w');

    % plot simulation results
    for j=1:nrRuns
        plot(simRes_1{j},projDimDyn{plotRun},'k');
        plot(simRes_2{j},projDimDyn{plotRun},'k');
        plot(simRes_3{j},projDimDyn{plotRun},'k');
    end


    % label plot
    xlabel(['x_{',num2str(projDimDyn{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projDimDyn{plotRun}(2)),'}']);
end

%% plot constraint variables
% projections
projDimAlg{1} = [1 31];
projDimAlg{2} = [2 32];
projDimAlg{3} = [3 33];
projDimAlg{4} = [4 34];
projDimAlg{5} = [5 35];
projDimAlg{6} = [6 36];
projDimAlg{7} = [7 37];      
projDimAlg{8} = [8 38];
projDimAlg{9} = [9 39];
projDimAlg{10} = [10 40];
projDimAlg{11} = [11 41];
projDimAlg{12} = [12 42];
projDimAlg{13} = [13 43];      
projDimAlg{14} = [14 44];   
projDimAlg{15} = [1 2];
projDimAlg{16} = [3 4];
projDimAlg{17} = [5 6];
projDimAlg{18} = [7 8];      

for plotRun=1:length(projDimAlg)
    
    h_alg{plotRun} = figure; hold on; box on;
    
    plotOrder = 10;
    
    %plot reachable sets of zonotope; mode 4
    for i=1:length(R_4.timeInterval.algebraic)
        Zproj = project(R_4.timeInterval.algebraic{i},projDimAlg{plotRun});
        Zproj = reduce(Zproj,'girard',plotOrder);
        plot(Zproj,[1 2],'FaceColor',[.6 .6 .6]);
    end
    
    %plot reachable sets of zonotope; mode 3
    for i=1:length(R_3.timeInterval.algebraic)
        Zproj = project(R_3.timeInterval.algebraic{i},projDimAlg{plotRun});
        Zproj = reduce(Zproj,'girard',plotOrder);
        plot(Zproj,[1 2],'FaceColor',[.6 .6 .6]);
    end 

    %plot reachable sets of zonotope; mode 1
    for i=1:length(R_1.timeInterval.algebraic)
        Zproj = project(R_1.timeInterval.algebraic{i},projDimAlg{plotRun});
        Zproj = reduce(Zproj,'girard',plotOrder);
        plot(Zproj,[1 2],'FaceColor',[.6 .6 .6]);
    end
    
    %plot reachable sets of zonotope; mode 2
    for i=1:length(R_2.timeInterval.algebraic)
        Zproj = project(R_2.timeInterval.algebraic{i},projDimAlg{plotRun});
        Zproj = reduce(Zproj,'girard',plotOrder);
        plot(Zproj,[1 2],'FaceColor',[.8 .8 .8]);
    end   
    
    % plot initial set
    plot(R_1.timeInterval.algebraic{1},projDimAlg{plotRun},'EdgeColor','k','FaceColor','w');

    % plot simulation results
    for j=1:nrRuns
        % mode 1
        for i=1:length(simRes_1{j}.t)
            plot(simRes_1{j}.a{1}(:,projDimAlg{plotRun}(1)), simRes_1{j}.a{1}(:,projDimAlg{plotRun}(2)),'Color',0*[1 1 1]);
        end
        % mode 2
        for i=1:length(simRes_2{j}.t)
            plot(simRes_2{j}.a{1}(:,projDimAlg{plotRun}(1)), simRes_2{j}.a{1}(:,projDimAlg{plotRun}(2)),'Color',0*[1 1 1]);
        end
        % mode 3
        for i=1:length(simRes_3{j}.t)
            plot(simRes_3{j}.a{1}(:,projDimAlg{plotRun}(1)), simRes_3{j}.a{1}(:,projDimAlg{plotRun}(2)),'Color',0*[1 1 1]);
        end
    end
    
    % label plot
    xlabel(['y_{',num2str(projDimAlg{plotRun}(1)),'}']);
    ylabel(['y_{',num2str(projDimAlg{plotRun}(2)),'}']);

end


%% save plots
cd(coraroot);
cd('../CORAplots/powerSystems/23Aug2022');

% dynamic variables
for iProj = 1:length(projDimDyn)
    %current projection
    currProj = projDimDyn{iProj};
    % file name
    fileName = ['IEEE30_dim_',num2str(currProj(1)),'_',num2str(currProj(2)),'_',datestr(now,'yyyymmmmdd_HH_MM_SS')];
    fileName = strrep(fileName,'.','_');
    
    % save
    saveas(h_dyn{iProj},[fileName,'.fig']);
end


% algebraic variables
for iProj = 1:length(projDimAlg)
    %current projection
    currProj = projDimAlg{iProj};
    % file name
    fileName = ['IEEE30alg_dim_',num2str(currProj(1)),'_',num2str(currProj(2)),'_',datestr(now,'yyyymmmmdd_HH_MM_SS')];
    fileName = strrep(fileName,'.','_');
    
    % save
    saveas(h_alg{iProj},[fileName,'.fig']);
end

% example completed
completed = 1;


% ------------------------------ END OF CODE ------------------------------
