function completed = example_nonlinear_reach_01_tank
% example_nonlinear_reach_01_tank - example of nonlinear reachability 
% analysis; this example is also a unit test function.
%
% This example can be found in Sec. 3.4.5 of
% M. Althoff, “Reachability analysis and its application to the safety 
% assessment of autonomous cars”, Dissertation, Technische Universität 
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
%    example_nonlinear_reach_01_tank
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
% Written:      18-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim=6;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=400; %final time
options.x0=[2; 4; 4; 2; 10; 4]; %initial state for simulation
options.R0=zonotope([options.x0,0.2*eye(dim)]); %initial state for reachability analysis

options.timeStep=1; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=50; %zonotope order
options.intermediateOrder=5;
options.reductionTechnique='girard';
options.errorOrder=1;
options.polytopeOrder=2; %polytope order
options.reductionInterval=1e3;
options.maxError = 1*ones(dim,1);

options.plotType='frame';
options.projectedDimensions=[1 2];

options.originContained = 0;
options.advancedLinErrorComp = 0;
options.tensorOrder = 2;
%--------------------------------------------------------------------------


%obtain uncertain inputs
options.uTrans = 0;
options.U = zonotope([0,0.005]); %input for reachability analysis

%specify continuous dynamics-----------------------------------------------
tank = nonlinearSys(6,1,@tank6Eq,options); %initialize tank system
%--------------------------------------------------------------------------


%compute reachable set using zonotopes
tic
Rcont = reach(tank, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(tank, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plot results--------------------------------------------------------------
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
    
    %plot reachable sets
    for i=1:length(Rcont)
        for j=1:length(Rcont{i})
            plotFilled(Rcont{i}{j},projectedDimensions,[.8 .8 .8],'EdgeColor','none');
        end
    end
    
    
    %plot initial set
    plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
    
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end

    %label plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------
