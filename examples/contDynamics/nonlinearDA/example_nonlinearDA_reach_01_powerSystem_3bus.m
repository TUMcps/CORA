function completed = example_nonlinearDA_reach_01_powerSystem_3bus()
% example_nonlinearDA_reach_01_powerSystem_3bus - example of 
% nonlinear-differntial-algebraic reachability analysis
%
% Computes the solution of the nonlinearDASys class for a 3 bus power system example;
%
% Syntax:  
%    example_nonlinearDA_reach_01_powerSystem_3bus()
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

options.tensorOrder = 1;

%specify continuous dynamics-----------------------------------------------
powerDyn = nonlinDASys(2,6,2,@bus3Dyn,@bus3Con,options); %initialize power system dynamics
%--------------------------------------------------------------------------

%set options --------------------------------------------------------------
options.tStart = 0; %start time
options.tFinal = 5; %final time
options.x0 = [380; 0.7]; %initial state for simulation (consistemcy of algebraic state is taken care of automatically)
options.y0guess = [ones(0.5*powerDyn.nrOfConstraints, 1); zeros(0.5*powerDyn.nrOfConstraints, 1)];
options.R0 = zonotope([options.x0,diag([0.1, 0.01])]); %initial state for reachability analysis
options.uTrans = [1; 0.4];
options.U = zonotope([zeros(2,1),diag([0, 0.1*options.uTrans(2)])]);

%options.timeStep=0.01; %time step size for reachable set computation
options.timeStep=0.05; %time step size for reachable set computation
options.taylorTerms=6; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order
options.errorOrder=1.5;
options.polytopeOrder=2; %polytope order
options.reductionTechnique='girard';

options.originContained = 0;
options.reductionInterval = 1e5;
options.advancedLinErrorComp = 0;

options.maxError = [0.5; 0];
options.maxError_x = options.maxError;
options.maxError_y = 0.005*[1; 1; 1; 1; 1; 1];
%--------------------------------------------------------------------------

%compute reachable set
tic
Rcont = reach(powerDyn, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(powerDyn, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plot results--------------------------------------------------------------
projectedDimensions=[1 2];
    
figure;
hold on

%plot reachable sets
for i=1:length(Rcont)
    for j=1:length(Rcont{1})
        Zproj = project(Rcont{i}{j},projectedDimensions);
        Zproj = reduce(Zproj,'girard',3);
        plotFilled(Zproj,[1 2],[.75 .75 .75],'EdgeColor','none');
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
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------
        