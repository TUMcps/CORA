function completed = example_hybrid_reach_ARCH18_drivetrain100p()
% example_hybrid_reach_02_powerTrain - example for hybrid dynamics; this 
% example is also a unit test function.
%
% This example can be found in
% M. Althoff and B. H. Krogh. Avoiding geometric intersection operations in 
% reachability analysis of hybrid systems. In Hybrid Systems: Computation 
% and Control, pages 45-54, 2012.
%
% Syntax:  
%    example_hybrid_reach_02_powerTrain()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      21-September-2011
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


% generate model
HA = drivetrain100p();

%set options---------------------------------------------------------------
Zcenter = [-0.0432;-11;0;30;0;30;360;-0.00132;30;0];
Zdelta = [0.0056;4.6666666666667;0;10;0;10;120;0.00056;10;0];
options.R0 = zonotope([Zcenter, Zdelta]); %initial state for reachability analysis
options.x0 = center(options.R0); %initial state for simulation
options.startLoc = 4; %initial location
options.finalLoc = 0; %0: no final location
options.timeStepLoc{1} = 5e-4; %time step size for reachable set computation
options.timeStepLoc{2} = 5e-4; %time step size for reachable set computation
options.timeStepLoc{3} = 5e-4; %time step size for reachable set computation
options.timeStepLoc{4} = 5e-4; %time step size for reachable set computation

options.taylorTerms = 20;
options.zonotopeOrder = 20;
options.errorOrder=1.2;
options.maxProjectionError=inf;

options.reductionTechnique = 'girard';
options.isHyperplaneMap = 0;
options.enclosureEnables = [3,5];
options.originContained = 0;
%--------------------------------------------------------------------------


%set u
u = [-5; 0]; %acceleration and load torque

%input:
for i = 1:4
    options.uLoc{i} = 0;
    options.uLocTrans{i} = 0;
    options.Uloc{i} = zonotope([0, 0]); 
end

%set start and end time
options.tStart = 0; %start time
options.tFinal = 2; %final time


%obtain random simulation results
for i=1:10
    %set initial state, input
    if i == 1
        options.x0 = Zcenter + Zdelta;
    elseif i == 2
        options.x0 = Zcenter - Zdelta; %initial state for simulation
    else
        options.x0 = randPoint(options.R0); %initial state for simulation
    end 
    
    %simulate hybrid automaton
    HAsim = simulate(HA,options); 
    simRes{i} = get(HAsim,'trajectory');
end

%compute reachable set
%with hyperplane intersection
tStart=tic;
%[HA_1] = reach_HSCC12(HA,options);
[HA_1] = reach(HA,options);
Rcont = get(HA_1,'continuousReachableSet');
tElapsed=toc(tStart)


for iPlot = 1:2
    
    %chhose projection
    if iPlot==1
        options.projectedDimensions = [1 2];
    elseif iPlot==2
        options.projectedDimensions = [1 3];
    end

    figure
    hold on
%     %plot continuous reachable set
%     redOrder = 3;
%     %1st time interval
%     Rplot = Rcont.OT;
%     %Rplot_noHP = Rcont_noHP.OT;
%     for iMode = 1:length(Rplot)
%         for iSet = 1:length(Rplot{iMode})
%             Rtmp = reduce(project(Rplot{iMode}{iSet},options.projectedDimensions),'girard',redOrder);
%             plotFilled(Rtmp,[1 2],[.6 .6 .6],'EdgeColor','none');            
%         end
%     end
%     
%     %2nd time interval
%     Rplot = Rcont2.OT;
%     %Rplot_noHP = Rcont_noHP2.OT;
%     for iMode = 1:length(Rplot)
%         for iSet = 1:length(Rplot{iMode})
%             Rtmp = reduce(project(Rplot{iMode}{iSet},options.projectedDimensions),'girard',redOrder);
%             plotFilled(Rtmp,[1 2],[.6 .6 .6],'EdgeColor','none');
%         end
%     end

    %plot initial set
    dim = length(options.x0);
    epsilonZ = zonotope(interval(-1e-4*ones(dim,1), 1e-4*ones(dim,1)));
    plotFilled(options.R0 + epsilonZ,options.projectedDimensions,'w','EdgeColor','k');

    %plot simulation results
    for n=1:length(simRes)
        for j=1:length(simRes{n}.x)
            %plot results
            plot(simRes{n}.x{j}(:,options.projectedDimensions(1)), simRes{n}.x{j}(:,options.projectedDimensions(2)),'k');
        end
    end
end

%example completed
completed = 1;

%------------- END OF CODE --------------