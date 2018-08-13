function completed = example_nonlinear_reach_05_autonomousCar()
% example_nonlinear_reach_05_autonomousCar - example of 
% nonlinear reachability analysis for following a reference trajectory; 
% this example is also a unit test function.
%
% This example is similar to the one in
%
% M. Althoff and J. M. Dolan. Online verification of automated road 
% vehicles using reachability analysis. IEEE Transactions on Robotics, 
% 30(4):903-918, 2014.
%
% One of the differences is that computational accelerations such as 
% parallelization are not activated for this example. Further
% accelerations, such as taking advantage of monotonicity in the Lagrange
% remainder are also not considered.
%
% Syntax:  
%    example_nonlinear_reach_05_autonomousCar
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


dim = 8;

options.maxError = ones(dim,1); % for comparison reasons


%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=3.99; %final time
options.x0=[0; 0; 0; 22; 0 ; 0; -2.1854; 0]; %initial state for simulation
options.R0 = zonotope([options.x0, 0.05*diag([1, 1, 1, 1, 1, 1, 1, 1])]); %initial state for reachability analysiszonotope([options.x0, diag([0.20, 0.20])]); %max for 3rd order
options.timeStep=0.01; %time step size for reachable set computation
options.taylorTerms=5; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order

options.uTransVec = uTRansVec4CASreach();
options.u = 0;
options.U = zonotope([0*options.uTransVec(:,1), 0.05*diag([ones(5,1);zeros(21,1)])]);

options.advancedLinErrorComp = 0;
options.tensorOrder = 1;
options.reductionInterval = inf;
options.reductionTechnique = 'girard';
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
vehicle = nonlinearSys(8,26,@vmodel_A_bicycle_linear_controlled,options); %initialize van-der-Pol oscillator
%--------------------------------------------------------------------------

%compute reachable set 
tic
Rcont = reach(vehicle, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(vehicle, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plot results--------------------------------------------------------------
for plotRun=1:3
    % plot different projections
    if plotRun==1
        projectedDimensions=[1 2];
        projectedReference=[17, 18];
    elseif plotRun==2
        projectedDimensions=[3 4];   
        projectedReference=[19, 20];
    elseif plotRun==3
        projectedDimensions=[5 6]; 
        projectedReference=[21, 22];
    end 
    
    figure;
    hold on

    %plot reachable sets 
    for i=1:length(Rcont)
        Zproj = project(Rcont{i}{1},projectedDimensions);
        Zproj = reduce(Zproj,'girard',3);
        plotFilled(Zproj,[1 2],[.8 .8 .8],'EdgeColor','none');
    end
    
    %plot initial set
    plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
    
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end
    
    %plot reference trajectory
    plot(options.uTransVec(projectedReference(1),:),options.uTransVec(projectedReference(2),:),'r','LineWidth',2);

    %label plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------

