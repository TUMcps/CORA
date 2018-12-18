function completed = example_nonlinearParam_reach_01_tank()
% example_nonlinearParam_reach_01_tank - example of nonlinear
% reachability analysis with uncertain parameters; this example is also a unit test function.
%
% This example can be found in Sec. 3.4.5 of
% M. Althoff, “Reachability analysis and its application to the safety 
% assessment of autonomous cars�?, Dissertation, Technische Universität 
% München, 2010, 
% http://nbn-resolving.de/urn/resolver.pl?urn:nbn:de:bvb:91-diss-20100715-963752-1-4
% 
% or in
%
% M. Althoff, O. Stursberg, and M. Buss. Reachability analysis of nonlinear 
% systems with uncertain parameters using conservative linearization. In 
% Proc. of the 47th IEEE Conference on Decision and Control, 
% pages 4042–4048, 2008
%
% Syntax:  
%    example_nonlinearParam_reach_01_tank
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
% Written:      19-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim=6;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=400; %final time
options.x0=[2; 4; 4; 2; 10; 4]; %initial state for simulation
options.R0=zonotope([options.x0,0.2*eye(dim)]); %initial state for reachability analysis
options.timeStep=1;
options.taylorTerms=4; %number of taylor terms for reachable sets
options.intermediateOrder = options.taylorTerms;
options.zonotopeOrder=10; %zonotope order
options.reductionTechnique='girard';
options.maxError = 1*ones(dim,1);
options.reductionInterval=1e3;
options.tensorOrder = 2;

options.advancedLinErrorComp = 0;

options.u=0; %input for simulation
options.U=zonotope([0,0.005]); %input for reachability analysis
options.uTrans=0; 

options.p=0.015; %parameter values for simulation
options.paramInt=interval(0.0148,0.015); %parameter intervals for reachability analysis
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%specify continuous dynamics with and without uncertain parameters---------
tankParam = nonlinParamSys(6,1,1,@tank6paramEq,options); %with uncertain parameters
tank = nonlinearSys(6,1,@tank6Eq,options); %without uncertain parameters
%--------------------------------------------------------------------------
        
%compute reachable set of tank system with and without uncertain parameters
tic
RcontParam = reach(tankParam,options); %with uncertain parameters
tComp = toc;
disp(['computation time of reachable set with uncertain parameters: ',num2str(tComp)]);
tic
RcontNoParam = reach(tank, options); %without uncertain parameters
tComp = toc;
disp(['computation time of reachable set without uncertain parameters: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(tank, options, runs, fractionVertices, fractionInputVertices, inputChanges);


%plot results--------------------------------------------------------------
plotOrder = 8;
for plotRun=1:3
    % plot different projections
    if plotRun==1
        projectedDimensions=[1 2];
    elseif plotRun==2
        projectedDimensions=[3 4];    
    elseif plotRun==3
        projectedDimensions=[5 6]; 
    end 

    figure;
    hold on

    %plot reachable sets of zonotope; uncertain parameters
    for i=1:length(RcontParam)
        for j=1:length(RcontParam{i})
            Zproj = reduce(RcontParam{i}{j},'girard',plotOrder);
            plotFilled(Zproj,projectedDimensions,[.675 .675 .675],'EdgeColor','none');
        end
    end
    
    %plot reachable sets of zonotope; without uncertain parameters
    for i=1:length(RcontNoParam)
        for j=1:length(RcontNoParam{i})
            Zproj = reduce(RcontNoParam{i}{j},'girard',plotOrder);
            plotFilled(Zproj,projectedDimensions,'w','EdgeColor','k');
        end
    end
    
    %plot initial set
    plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
    
  
    %plot simulation results      
    for i=1:length(simRes.x)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'k');
    end

    %label plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------
