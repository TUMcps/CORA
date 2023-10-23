function res = testLong_nonlinearSys_reach_05_autonomousCar()
% testLong_nonlinearSys_reach_05_autonomousCar - unit_test_function of 
%    nonlinear reachability analysis for following a reference trajectory
%    Checks the solution of an autonomous car following a reference
%    trajectory; It is checked whether the reachable set is enclosed
%    in the initial set after a certain amount of time.
%
% Syntax:
%    res = testLong_nonlinearSys_reach_05_autonomousCar()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       10-September-2015
% Last update:   12-August-2016
%                23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

dim_x = 8;
params.tFinal=0.1; %final time
params.R0 = zonotope([[0; 0; 0; 22; 0 ; 0; -2.1854; 0],...
    0.05*diag([1, 1, 1, 1, 1, 1, 1, 1])]); %initial state for reachability analysis
params.u = uTRansVec4CASreach();
params.u = params.u(:,1:10);
params.U = zonotope([0*params.u(:,1), 0.05*diag([ones(5,1);zeros(21,1)])]);


% Reachability Settings ---------------------------------------------------

options.timeStep=0.01; %time step size for reachable set computation
options.taylorTerms=5; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order
options.maxError = ones(dim_x,1); % for comparison reasons

options.alg = 'lin';
options.tensorOrder = 2;
options.reductionInterval = inf;


% System Dynamics ---------------------------------------------------------

vehicle = nonlinearSys(@vmodel_A_bicycle_linear_controlled);


% Reachability Analysis --------------------------------------------------- 

R = reach(vehicle, params, options);


% Numerical Evaluation ----------------------------------------------------

%enclose result by interval
IH = interval(R.timeInterval.set{end});

%saved result
IH_saved = interval( ...
    [1.9113165970926538; -0.1763919452055745; -0.0628054832382009; 21.7144766634198660; -0.1081824122890472; -0.2077292973551603; -2.5375737849902964; -0.0418129344907171], ...
    [2.2486917068056300; 0.1775907109675143; 0.0638901159110268; 21.8628568504537739; 0.1329959480608928; 0.2499135775758194; -1.9819546806525596; 0.0646318796537233]);

%final result
res = isequal(IH,IH_saved,1e-8);

% %simulate
% stepsizeOptions = odeset('MaxStep',0.2*(options.tFinal-options.tStart));
% %generate overall options
% opt = odeset(stepsizeOptions);
% 
% 
% %initialize
% runs=40;
% inputChanges=length(Rcont);
% t=cell(runs,inputChanges);
% x=cell(runs,inputChanges);
% finalTime=options.tFinal;
% R0 = reduce(options.R0,'girard',1);
% 
% for i=1:runs
%     params.tStart=0;
%     params.tFinal=finalTime/inputChanges;
%     for iChange = 1:inputChanges
%         %set initial state, input
%         if iChange == 1
%             if i<=30
%                 params.x0=randPoint(R0,1,'extreme'); %initial state for simulation
%             else
%                 params.x0=randPoint(R0); %initial state for simulation
%             end
%         else
%             params.tStart=params.tFinal;
%             params.tFinal=params.tFinal+finalTime/inputChanges;
%             params.x0 = x{i,iChange-1}(end,:);
%         end
%         
%         %set input
%         params.uTrans = params.uTransVec(:,iChange);
%         if i<=8
%             params.u=randPoint(options.U,1,'extreme')+params.uTrans; %input for simulation
%         else
%             params.u=randPoint(options.U)+params.uTrans; %input for simulation
%         end
% 
%         %simulate hybrid automaton
%         [vehicle,t{i,iChange},x{i,iChange}] = simulate(vehicle,params,opt); 
%     end
% end
% 
% 
% plotOrder = 10;
% 
% %plot dynamic variables
% for plotRun=1:2
% 
%     if plotRun==1
%         projectedDimensions=[1 2];
%     elseif plotRun==2
%         projectedDimensions=[3 4];     
%     end 
%     
% 
%     figure;
%     hold on
% 
%     %plot reachable sets of zonotope; mode 1
%     for i=1:length(Rcont)
%         Zproj = project(Rcont{i}{1},projectedDimensions);
%         Zproj = reduce(Zproj,'girard',plotOrder);
%         plot(Zproj,[1 2],'FaceColor',[.75 .75 .75],'Filled',true,'EdgeColor','none');
%     end
%     
% %     %plot reachable sets of zonotope; mode 1
% %     for i=2:length(Rcont)
% %         plot(Rcont{i},projectedDimensions,'lightgray',2,plotOrder);
% %     end
%     
% %     pSet = randPoint(project(Rcont{end},[1 2]), 1e3);
% %     plot(pSet(1,:), pSet(2,:), 'r.');
%     
%     %plot simulation results      
%     for i=1:length(t(:,1))
%         for j=1:length(t(i,:))
%             plot(x{i,j}(:,projectedDimensions(1)),x{i,j}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
%         end
%     end
%     
%     %plot initial set
%     plot(R0,projectedDimensions,'w-','lineWidth',2);
% 
%     xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
%     ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
% end

% ------------------------------ END OF CODE ------------------------------
