function completed = example_hybrid_reach_ARCH18_powerTrain_DTN04()
% example_hybrid_reach_ARCH18_powerTrain_DTN04 - example for hybrid 
% dynamics from the ARCH18 friendly competition (Test case DTN04)
%
% This example can be found in
% M. Althoff and B. H. Krogh. Avoiding geometric intersection operations in 
% reachability analysis of hybrid systems. In Hybrid Systems: Computation 
% and Control, pages 45-54, 2012.
%
% Syntax:  
%    example_hybrid_reach_ARCH18_powerTrain_DTN04()
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
% Author:       Matthias Althof, Niklas Kochdumper
% Written:      21-September-2011
% Last update:  29-May-2018
% Last revision:---


%------------- BEGIN CODE --------------


% generate model
[HA, Zcenter, Zdelta, B1, B2, B3, c1, c2, c3] = initPowerTrain(51); 
Zdelta = 0.05*Zdelta;



% Options -----------------------------------------------------------------

options.R0 = zonotope([Zcenter, Zdelta]); %initial state for reachability analysis
options.x0 = center(options.R0); %initial state for simulation
options.startLoc = 3; %initial location
options.finalLoc = 0; %0: no final location
options.timeStepLoc{1} = 5e-4; %time step size for reachable set computation
options.timeStepLoc{2} = 5e-4; %time step size for reachable set computation
options.timeStepLoc{3} = 5e-4; %time step size for reachable set computation

options.taylorTerms = 20;
options.zonotopeOrder = 20;
options.errorOrder=1.2;
options.maxProjectionError=inf;
options.linAlg = 1;

options.reductionTechnique = 'girard';
options.isHyperplaneMap = 1;
options.originContained = 0;
options.guardIntersect = 'hyperplaneMap';


% System input
u = [-5; 0]; %acceleration and load torque

%input:
%loc 1
options.uLoc{1} = B1*u + c1;
options.uLocTrans{1} = options.uLoc{1};
options.Uloc{1} = B1*zonotope([0; 0]); 
%loc 2
options.uLoc{2} = B2*u + c2;
options.uLocTrans{2} = options.uLoc{2};
options.Uloc{2} = B2*zonotope([0; 0]); 
%loc 3
options.uLoc{3} = B3*u + c3;
options.uLocTrans{3} = options.uLoc{3};
options.Uloc{3} = B3*zonotope([0; 0]); 

%set start and end time
options.tStart = 0; %start time
options.tFinal = 0.2; %final time




% Simulation 11------------------------------------------------------------

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




% Reachability Analysis 1 -------------------------------------------------

tStart=tic;
%[HA_1] = reach_HSCC12(HA,options);
[HA_1] = reach(HA,options);
Rcont = get(HA_1,'continuousReachableSet');
disp('FIRST PART')
tElapsed=toc(tStart)



% Options -----------------------------------------------------------------

options.tStart = 0.2; %start time
options.tFinal = 2; %final time
%set u
u = [5; 0]; %acceleration and load torque

%input:
%loc 1
options.uLoc{1} = B1*u + c1;
options.uLocTrans{1} = options.uLoc{1};
options.Uloc{1} = B1*zonotope([0; 0]); 
%loc 2
options.uLoc{2} = B2*u + c2;
options.uLocTrans{2} = options.uLoc{2};
options.Uloc{2} = B2*zonotope([0; 0]); 
%loc 3
options.uLoc{3} = B3*u + c3;
options.uLocTrans{3} = options.uLoc{3};
options.Uloc{3} = B3*zonotope([0; 0]); 



% Simulation 2 ------------------------------------------------------------

%obtain random simulation results
for i=1:10
    %set initial state
    options.x0=simRes{i}.x{1}(end,:); %initial state for simulation
    
    %simulate hybrid automaton
    HAsim = simulate(HA,options); 
    simRes2{i} = get(HAsim,'trajectory');
end



% Reachability Analysis 2 -------------------------------------------------

%profile on
tStart=tic;
options.R0 = Rcont.T{end}{end};
%[HA_2] = reach_HSCC12(HA,options);
[HA_2] = reach(HA,options);
Rcont2 = get(HA_2,'continuousReachableSet');
Rencl = get(HA_2,'enclosure');
disp('SECOND PART')
tElapsed=toc(tStart)




% Visualization -----------------------------------------------------------

for iPlot = 1:2
    
    %chhose projection
    if iPlot==1
        options.projectedDimensions = [1 2];
        xlab = 'x_1';
        ylab = 'x_2';
    elseif iPlot==2
        options.projectedDimensions = [1 3];
        xlab = 'x_1';
        ylab = 'x_3';
    end

    figure
    hold on
    %plot continuous reachable set
    redOrder = 3;
    %1st time interval
    Rplot = Rcont.OT;
    %Rplot_noHP = Rcont_noHP.OT;
    for iMode = 1:length(Rplot)
        for iSet = 1:length(Rplot{iMode})
            Rtmp = reduce(project(Rplot{iMode}{iSet},options.projectedDimensions),'girard',redOrder);
            plotFilled(Rtmp,[1 2],[.6 .6 .6],'EdgeColor','none');            
        end
    end
    
    %2nd time interval
    Rplot = Rcont2.OT;
    %Rplot_noHP = Rcont_noHP2.OT;
    for iMode = 1:length(Rplot)
        for iSet = 1:length(Rplot{iMode})
            Rtmp = reduce(project(Rplot{iMode}{iSet},options.projectedDimensions),'girard',redOrder);
            plotFilled(Rtmp,[1 2],[.6 .6 .6],'EdgeColor','none');
        end
    end

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
    for n=1:length(simRes2)
        for j=1:length(simRes2{n}.x)
            %plot results
            plot(simRes2{n}.x{j}(:,options.projectedDimensions(1)), simRes2{n}.x{j}(:,options.projectedDimensions(2)),'k');
        end
    end
    
    box on
    xlabel(xlab);
    ylabel(ylab);

end

%example completed
completed = 1;

%------------- END OF CODE --------------
