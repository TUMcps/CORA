function completed = example_nonlinearDA_reach_05_powerSystem_14bus_compositional()
% example_nonlinearDA_reach_05_powerSystem_14bus_compositional - example of 
%    nonlinear-differential-algebraic reachability analysis: 14-bus power 
%    system from [1], where the abstraction error is computed 
%    compositionally.
%
% Syntax:
%    completed = example_nonlinearDA_reach_05_powerSystem_14bus_compositional()
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
% Written:       08-June-2022
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
    load([path filesep 'IEEE14_model'], 'IEEE14_model');
    derivatives(IEEE14_model, options); % derivatives are not automatically created
end
if ~isfile([path filesep 'IEEE14_sub1_model.mat'])
    powerSystem2cora('IEEE14_sub1');
    load([path filesep 'IEEE14_sub1_model'], 'IEEE14_sub1_model');
    derivatives(IEEE14_sub1_model, options); % derivatives are not automatically created for subsystems
end
if ~isfile([path filesep 'IEEE14_sub2_model.mat'])
    powerSystem2cora('IEEE14_sub2');
    load([path filesep 'IEEE14_sub2_model'], 'IEEE14_sub2_model');
    derivatives(IEEE14_sub2_model, options); % derivatives are not automatically created for subsystems
end
if ~isfile([path filesep 'IEEE14_fault_model.mat'])
    powerSystem2cora('IEEE14_fault');
    load([path filesep 'IEEE14_fault_model'], 'IEEE14_fault_model');
    derivatives(IEEE14_model, options); % derivatives are not automatically created
end
if ~isfile([path filesep 'IEEE14_fault_sub1_model.mat'])
    powerSystem2cora('IEEE14_fault_sub1');
    load([path filesep 'IEEE14_fault_sub1_model'], 'IEEE14_fault_sub1_model');
    derivatives(IEEE14_fault_sub1_model, options); % derivatives are not automatically created for subsystems
end

% load models (normal and faulty operation)
load([path filesep 'IEEE14_model'], 'IEEE14_model');
load([path filesep 'IEEE14_sub1_model'], 'IEEE14_sub1_model');
load([path filesep 'IEEE14_sub2_model'], 'IEEE14_sub2_model');
load([path filesep 'IEEE14_fault_model'], 'IEEE14_fault_model');
load([path filesep 'IEEE14_fault_sub1_model'], 'IEEE14_fault_sub1_model');

% rewrite systems as cells
sys{1} = IEEE14_model;
sys{2} = IEEE14_sub1_model;
sys{3} = IEEE14_sub2_model;
sys_fault{1} = IEEE14_fault_model;
sys_fault{2} = IEEE14_fault_sub1_model;
sys_fault{3} = IEEE14_sub2_model; % this subsystem has no fault

% build indices for the input, state, and constraints relating the full 
% full system
load IEEE14.mat bus
busFull = bus;
% subsystem 1
load IEEE14_sub1.mat bus
bus_subsystem{1} = bus;
% subsystem 2
load IEEE14_sub2.mat bus
bus_subsystem{2} = bus;
% create cell of subsystems
subsystem{1} = sys{2};
subsystem{2} = sys{3};
% compute indices for subsystems
index = indexForSubsystems(sys{1}, subsystem, busFull, bus_subsystem, 3);

% faulty system: build indices for the input, state, and constraints 
% relating the full 
% full system
load IEEE14_fault.mat bus
busFull = bus;
% subsystem 1
load IEEE14_fault_sub1.mat bus
bus_subsystem{1} = bus;
% subsystem 2; this subsystem has no fault
load IEEE14_sub2.mat bus
bus_subsystem{2} = bus;
% create cell of subsystems
subsystem{1} = sys_fault{2};
subsystem{2} = sys_fault{3};
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
load IEEE14.mat VM bus
options.Vgen = VM(bus.generator);
% store Jacobian and Hessians of subsystems
sys{2}.setHessian('int');
sys{3}.setHessian('int');
options.subsystemHessian{1} = sys{2}.hessian;
options.subsystemHessian{2} = sys{3}.hessian;
options.subsystemJacobian{1} = sys{2}.jacobian;
options.subsystemJacobian{2} = sys{3}.jacobian;
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
% provide index for compositional computation
options.index = index_fault;
% load model parameters to obtain Vgen
load IEEE14_fault.mat VM bus
options.Vgen = VM(bus.generator);
% store Jacobian and Hessians of subsystems
sys_fault{2}.setHessian('int');
sys_fault{3}.setHessian('int');
options.subsystemHessian{1} = sys_fault{2}.hessian;
options.subsystemHessian{2} = sys_fault{3}.hessian;
options.subsystemJacobian{1} = sys_fault{2}.jacobian;
options.subsystemJacobian{2} = sys_fault{3}.jacobian;
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
timeSteps = round((params.tFinal - params.tStart)/options.timeStep);
% provide index for compositional computation
options.index = index;
% load model parameters to obtain Vgen
load IEEE14.mat VM bus
options.Vgen = VM(bus.generator);
% store Jacobian and Hessians of subsystems
options.subsystemHessian{1} = sys{2}.hessian;
options.subsystemHessian{2} = sys{3}.hessian;
options.subsystemJacobian{1} = sys{2}.jacobian;
options.subsystemJacobian{2} = sys{3}.jacobian;
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

%% operation after fixing the fault (larger time step size)
% time horizon and time step
params.tStart = 2; 
params.tFinal = 5; 
options.timeStep = 0.02; 
timeSteps = round((params.tFinal - params.tStart)/options.timeStep);
% input
u0_gen = [2; 0.4; 0; 0; 0];
params.u = [u0_gen; u0_load]*ones(1,timeSteps);
% initial set is final set of normal operation
params.R0 = R_3.timePoint.set{end}; %initial state for reachability analysis
% computation of reachable set
tic
R_4 = reach(IEEE14_model, params, options);
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
    params.tFinal = 5; 
    u0_gen = [2; 0.4; 0; 0; 0];
    params.u = [u0_gen; u0_load];
    % initial set is final solution of normal operation
    params.R0 = zonotope(simRes_2{iRun}.x{1}(end,:)');
    simRes_3{iRun} = simulateRandom(IEEE14_model, params, simOpt);
end
    
% Visualization -----------------------------------------------------------

%plot dynamic variables
for plotRun=1:19
    if plotRun==1
        projDim=[1 6];
    elseif plotRun==2
        projDim=[2 7];
    elseif plotRun==3
        projDim=[3 8];
    elseif plotRun==4
        projDim=[4 9];
    elseif plotRun==5
        projDim=[5 10];
    elseif plotRun==6
        projDim=[11 12];
    elseif plotRun==7
        projDim=[13 14];      
    elseif plotRun==8
        projDim=[6 11];
    elseif plotRun==9
        projDim=[7 12];
    elseif plotRun==10
        projDim=[8 13];
    elseif plotRun==11
        projDim=[9 14];      
    elseif plotRun==12
        projDim=[1 11];
    elseif plotRun==13
        projDim=[2 12];
    elseif plotRun==14
        projDim=[3 13];      
    elseif plotRun==15
        projDim=[4 14];   
    elseif plotRun==16
        projDim=[1 2];
    elseif plotRun==17
        projDim=[3 4];
    elseif plotRun==18
        projDim=[5 6];
    elseif plotRun==19
        projDim=[7 8];      
    end 

    figure; hold on; box on;

    % plot reachable sets
    plot(R_4, projDim, 'Order', 10, 'FaceColor',[.6 .6 .6]);
    plot(R_3, projDim, 'Order', 10, 'FaceColor',[.6 .6 .6]);
    plot(R_2, projDim, 'Order', 10, 'FaceColor',[.8 .8 .8]);
    plot(R_1, projDim, 'Order', 10, 'FaceColor',[.6 .6 .6]);
    
    
    % plot initial set
    plot(params.R0,projDim,'EdgeColor','k','FaceColor','w');


    % plot simulation results
    for j=1:nrRuns
        plot(simRes_1{j},projDim,'k');
        plot(simRes_2{j},projDim,'k');
        plot(simRes_3{j},projDim,'k');
    end


    % label plot
    xlabel(['x_{',num2str(projDim(1)),'}']);
    ylabel(['x_{',num2str(projDim(2)),'}']);
end

%plot constraint variables
for plotRun=1:18
    if plotRun==1
        projDim=[1 15];
    elseif plotRun==2
        projDim=[2 16];
    elseif plotRun==3
        projDim=[3 17];
    elseif plotRun==4
        projDim=[4 18];
    elseif plotRun==5
        projDim=[5 19];
    elseif plotRun==6
        projDim=[6 20];
    elseif plotRun==7
        projDim=[7 21];      
    elseif plotRun==8
        projDim=[8 22];
    elseif plotRun==9
        projDim=[9 23];
    elseif plotRun==10
        projDim=[10 24];
    elseif plotRun==11
        projDim=[11 25];
    elseif plotRun==12
        projDim=[12 26];
    elseif plotRun==13
        projDim=[13 27];      
    elseif plotRun==14
        projDim=[14 27];   
    elseif plotRun==15
        projDim=[1 2];
    elseif plotRun==16
        projDim=[3 4];
    elseif plotRun==17
        projDim=[5 6];
    elseif plotRun==18
        projDim=[7 8];      
    end 
    
    figure; hold on; box on;
    
    plotOrder = 10;
    
    %plot reachable sets of zonotope; mode 4
    for i=1:length(R_4.timeInterval.algebraic)
        Zproj = project(R_4.timeInterval.algebraic{i},projDim);
        Zproj = reduce(Zproj,'girard',plotOrder);
        plot(Zproj,[1 2],'FaceColor',[.6 .6 .6]);
    end
    
    %plot reachable sets of zonotope; mode 3
    for i=1:length(R_3.timeInterval.algebraic)
        Zproj = project(R_3.timeInterval.algebraic{i},projDim);
        Zproj = reduce(Zproj,'girard',plotOrder);
        plot(Zproj,[1 2],'FaceColor',[.6 .6 .6]);
    end 

    %plot reachable sets of zonotope; mode 1
    for i=1:length(R_1.timeInterval.algebraic)
        Zproj = project(R_1.timeInterval.algebraic{i},projDim);
        Zproj = reduce(Zproj,'girard',plotOrder);
        plot(Zproj,[1 2],'FaceColor',[.6 .6 .6]);
    end
    
    %plot reachable sets of zonotope; mode 2
    for i=1:length(R_2.timeInterval.algebraic)
        Zproj = project(R_2.timeInterval.algebraic{i},projDim);
        Zproj = reduce(Zproj,'girard',plotOrder);
        plot(Zproj,[1 2],'FaceColor',[.8 .8 .8]);
    end   

    % plot simulation results
    for j=1:nrRuns
        % mode 1
        for i=1:length(simRes_1{j}.t)
            plot(simRes_1{j}.a{1}(:,projDim(1)), simRes_1{j}.a{1}(:,projDim(2)),'Color',0*[1 1 1]);
        end
        % mode 2
        for i=1:length(simRes_2{j}.t)
            plot(simRes_2{j}.a{1}(:,projDim(1)), simRes_2{j}.a{1}(:,projDim(2)),'Color',0*[1 1 1]);
        end
        % mode 3
        for i=1:length(simRes_3{j}.t)
            plot(simRes_3{j}.a{1}(:,projDim(1)), simRes_3{j}.a{1}(:,projDim(2)),'Color',0*[1 1 1]);
        end
    end

    xlabel(['y_{',num2str(projDim(1)),'}']);
    ylabel(['y_{',num2str(projDim(2)),'}']);
    
    % label plot
    xlabel(['y_{',num2str(projDim(1)),'}']);
    ylabel(['y_{',num2str(projDim(2)),'}']);

end

% example completed
completed = 1;


% ------------------------------ END OF CODE ------------------------------
