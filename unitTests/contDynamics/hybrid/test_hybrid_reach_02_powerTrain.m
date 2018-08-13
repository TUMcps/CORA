function res = test_hybrid_reach_02_powerTrain()
% test_hybrid_reach_02_powerTrain - unit_test_function for hybrid
% dynamics
%
% Checks the solution of the hybrid system class for the powertrain example
% in HSCC'12;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been 
% saved. It is also checked whether the simulation matches the analytical
% solution.
%
% Syntax:  
%    res = test_hybrid_reach_02_powerTrain
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
[HA, Zcenter, Zdelta, B1, B2, B3, c1, c2, c3] = initPowerTrain(11); 

%set options---------------------------------------------------------------
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

options.reductionTechnique = 'girard';
options.isHyperplaneMap = 1;
options.originContained = 0;
options.linAlg = 1;
%--------------------------------------------------------------------------


%set u
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
disp('FIRST PART')
tElapsed=toc(tStart)

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

%obtain random simulation results
for i=1:10
    %set initial state
    options.x0=simRes{i}.x{1}(end,:); %initial state for simulation
    
    %simulate hybrid automaton
    HAsim = simulate(HA,options); 
    simRes2{i} = get(HAsim,'trajectory');
end

%compute reachable set
%with hyperplane intersection
%profile on
tStart=tic;
options.R0 = Rcont.T{end}{end};
%[HA_2] = reach_HSCC12(HA,options);
[HA_2] = reach(HA,options);
Rcont2 = get(HA_2,'continuousReachableSet');
Rencl = get(HA_2,'enclosure');
disp('SECOND PART')
tElapsed=toc(tStart)


%CHECK REACHABLE SET-------------------------------------------------------
R = get(HA_2,'continuousReachableSet');
I = interval(R.OT{end}{end});

%saved result
I_saved = interval( ...
[0.1182651639334430;83.0046825479405328;46.1583855538025603;28.0182406232294383;40.5392194787527842;21.5475780200804721;260.1116869253353912;40.5494659581651860;21.6287962190078389;40.5420330233180408;21.5815099154996162], ...
[0.1307501306623376;91.4773598749626018;87.6407011390990078;48.0963152143391355;81.3400364607803823;41.1397219665579357;494.6158000089853317;81.3529246583789813;41.2569184405087555;81.3470509834557021;41.3218975384758593]);
        
%check if slightly bloated versions enclose each other
res_1 = (I <= enlarge(I_saved,1+1e-8));
res_2 = (I_saved <= enlarge(I,1+1e-8));

%final result
res = res_1*res_2;
%--------------------------------------------------------------------------


% for iPlot = 1:2
%     
%     %chhose projection
%     if iPlot==1
%         options.projectedDimensions = [1 2];
%     elseif iPlot==2
%         options.projectedDimensions = [1 3];
%     end
% 
%     figure
%     hold on
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
% 
%     %plot initial set
%     dim = length(options.x0);
%     epsilonZ = zonotope(interval(-1e-4*ones(dim,1), 1e-4*ones(dim,1)));
%     plotFilled(options.R0 + epsilonZ,options.projectedDimensions,'w','EdgeColor','k');
% 
%     %plot simulation results
%     for n=1:length(simRes)
%         for j=1:length(simRes{n}.x)
%             %plot results
%             plot(simRes{n}.x{j}(:,options.projectedDimensions(1)), simRes{n}.x{j}(:,options.projectedDimensions(2)),'k');
%         end
%     end
%     for n=1:length(simRes2)
%         for j=1:length(simRes2{n}.x)
%             %plot results
%             plot(simRes2{n}.x{j}(:,options.projectedDimensions(1)), simRes2{n}.x{j}(:,options.projectedDimensions(2)),'k');
%         end
%     end
% 
% end

%------------- END OF CODE --------------
