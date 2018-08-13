function example_nonlinear_reach_ARCH18_quadrocopterControlled
% example_nonlinear_reach_ARCH18_quadrocopter - example of nonlinear reachability 
% analysis in the ARCH'18 friendly competition
%
% Syntax:  
%    example_nonlinear_reach_ARCH18_quadrocopterControlled
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
% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      15-April-2017
% Last update:  30-May-2018
% Last revision:---


%------------- BEGIN CODE --------------

dim=12;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=5; %final time
options.x0=zeros(dim,1); %initial state for simulation
options.R0=zonotope([options.x0,0.4*diag([1; 1; 1; 1; 1; 1; zeros(6,1)])]); %initial state for reachability analysis

%options.timeStep=0.025; %time step size for reachable set computation
%options.timeStep=0.05; %time step size for reachable set computation
options.timeStep=0.1; %time step size for reachable set computation
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
options.uTrans = [1;0;0];
options.U = zonotope([0;0;0]); %input for reachability analysis

%specify continuous dynamics-----------------------------------------------
quadrocopter = nonlinearSys(12,3,@quadrocopterControlledEq,options); %initialize tank system
%--------------------------------------------------------------------------


%compute reachable set using zonotopes
tic
Rcont = reach(quadrocopter, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
%runs = 3;
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(quadrocopter, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plot results--------------------------------------------------------------
for plotRun=1:1
    % plot different projections
    if plotRun==1
        projectedDimensions=[1 2];
    elseif plotRun==2
        projectedDimensions=[3 4];    
    elseif plotRun==3
        projectedDimensions=[5 6]; 
    elseif plotRun==4
        projectedDimensions=[7 8];
    end 
    
    figure;
    hold on

%     %plot reachable sets 
%     for i=1:length(Rcont)
%         plotFilled(Rcont{i}{1},projectedDimensions,[.8 .8 .8],'EdgeColor','none');
%     end

    %plot time elapse
    for i=1:length(Rcont)
        %get Uout 
        Uout1 = interval(project(Rcont{i}{1},3));
        %obtain times
        t1 = (i-1)*options.timeStep;
        t2 = i*options.timeStep;
        %generate plot areas as interval hulls
        IH = interval([t1; infimum(Uout1)], [t2; supremum(Uout1)]);

        plotFilled(IH,[1 2],[.75 .75 .75],'EdgeColor','none');
    end

    %plot simulation results
    for i=1:(length(simRes.t))
        plot(simRes.t{i},simRes.x{i}(:,3),'Color',0*[1 1 1]);
    end
    
    %label plot
    box on
    xlabel(['t']);
    ylabel(['altitude x_3']);
    axis([0,5,-0.6,1.5]);
    
    
%     %plot initial set
%     plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
%     
%     %plot simulation results      
%     for i=1:length(simRes.t)
%         plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
%     end
% 
%     %label plot
%     xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
%     ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end
%--------------------------------------------------------------------------

%------------- END OF CODE --------------
