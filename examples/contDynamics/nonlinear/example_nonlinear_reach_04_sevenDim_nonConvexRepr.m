function completed = example_nonlinear_reach_04_sevenDim_nonConvexRepr()
% example_nonlinear_reach_04_sevenDim_nonConvexRepr - example of 
% nonlinear reachability analysis; this example is also a unit test function.
%
% This example originates from an email exchange with Xin Chen; check if
% part of some benchmark suite.
%
% Syntax:  
%    example_nonlinear_reach_04_sevenDim_nonConvexRepr
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


dim=7;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=5; %final time
options.x0=[1.2; 1.05; 1.5; 2.4; 1; 0.1; 0.45]; %initial state for simulation
options.R0=quadZonotope(options.x0,0.2*eye(dim),[],[],[]); %initial state for reachability analysis

options.timeStep=0.005; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order
options.intermediateOrder=5;
options.reductionTechnique='girard';
options.errorOrder=2;
options.polytopeOrder=10; %polytope order
options.reductionInterval=inf;
options.maxError = 2*ones(dim,1);

options.plotType='frame';
options.projectedDimensions=[1 2];

options.originContained = 0;
options.advancedLinErrorComp = 1;
options.tensorOrder = 3;
%--------------------------------------------------------------------------


%obtain uncertain inputs
options.uTrans = 0;
options.U = zonotope([0]); %input for reachability analysis

%specify continuous dynamics-----------------------------------------------
sys=nonlinearSys(7,1,@sevenDimNonlinEq,options); %initialize system
%--------------------------------------------------------------------------


%compute reachable set using polynomial zonotopes
tic
Rcont = reach(sys, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(sys, options, runs, fractionVertices, fractionInputVertices, inputChanges);

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
             plotFilled(Rcont{i}{j},projectedDimensions,1,8,[.8 .8 .8],'EdgeColor','none');
        end
    end

    %plot initial set
    plotFilled(options.R0,projectedDimensions,1,[],'w','EdgeColor','k');

    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end

    %öabel plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------