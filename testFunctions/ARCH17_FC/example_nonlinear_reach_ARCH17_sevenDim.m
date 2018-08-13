function example_nonlinear_reach_ARCH17_sevenDim()
% example_nonlinear_reach_ARCH17_sevenDim - example of 
% nonlinear reachability analysis
%
% Syntax:  
%    example_nonlinear_reach_ARCH17_sevenDim
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
% Written:      26-March-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


dim=7;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=20; %final time
options.x0=[1.2; 1.05; 1.5; 2.4; 1; 0.1; 0.45]; %initial state for simulation
options.R0=quadZonotope(options.x0,0.1*eye(dim),[],[],[]); %initial state for reachability analysis

%options.timeStep=0.01; %time step size for reachable set computation
options.timeStep=0.05; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=50; %zonotope order
options.intermediateOrder=5;
options.reductionTechnique='girard';
options.errorOrder=2;
options.polytopeOrder=10; %polytope order
options.reductionInterval=1e3;
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

%plot results over time
for iDim = 4:4

    figure;
    hold on
    
    t1 = 0;
    t2 = options.timeStep;
    init = 1;

    %plot time elapse
    for i=1:length(Rcont)
        for j=1:length(Rcont{i}) 
            %get Uout 
            Uout1 = interval(project(Rcont{i}{j},iDim));
            %obtain times
            if init
                init = 0;
            else
                t1 = t1 + options.timeStep;
                t2 = t2 + options.timeStep;
            end
            %generate plot areas as interval hulls
            IH1 = interval([t1; infimum(Uout1)], [t2; supremum(Uout1)]);

            plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
        end
    end

    %plot simulation results
    for i=1:(length(simRes.t))
        plot(simRes.t{i},simRes.x{i}(:,iDim),'Color',0*[1 1 1]);
    end
end

xlabel('t')
ylabel('x_4')
axis([0,20,1.3,4.8])


%------------- END OF CODE --------------