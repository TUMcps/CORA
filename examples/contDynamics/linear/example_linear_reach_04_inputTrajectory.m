function completed = example_linear_reach_04_inputTrajectory()
% example_linear_reach_04_inputTrajectory - example for linear reachability 
% analysis with an input trajectory uTransVec; this example is also a unit
% test function.
%
% Syntax:  
%    completed = example_linear_reach_04_inputTrajectory()
%
% Inputs:
%    no
%
% Outputs:
%    completed - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      27-July-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% load data
load data_Jean-MarcBiannic.mat A B uvec
dim = length(A);
inputDim = length(B(1,:));

%set options --------------------------------------------------------------
options.tStart = 0; %start time
options.tFinal = 10; %final time
options.x0 = zeros(dim,1); %initial state for simulation
options.R0 = zonotope([options.x0,0.1*eye(dim,4)]); %initial state for reachability analysis

options.timeStep = 0.01; %time step size for reachable set computation
options.taylorTerms = 4; %number of taylor terms for reachable sets
options.zonotopeOrder = 50; %zonotope order
options.originContained = 0;
options.reductionTechnique = 'girard';
options.linAlg = 1;

options.uTransVec=[zeros(inputDim,1) uvec]; % input trajectory
options.U = zonotope([zeros(inputDim,1),diag([0.05 1])]); %input for reachability analysis

options.linAlg=2;
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
sys = linearSys('JeanMarcSys',A,B); %initialize system
%--------------------------------------------------------------------------

%compute reachable set using zonotopes
tic
Rcont = reach(sys, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 20;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 10;
simRes = simulate_random(sys, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plot results--------------------------------------------------------------
for plotRun=1:2
    % plot different projections
    if plotRun==1
        projectedDimensions=[1 2];
    elseif plotRun==2
        projectedDimensions=[3 4];     
    end 
    
    figure;
    hold on

    %plot reachable sets 
    for i=1:length(Rcont)
        Zproj = project(Rcont{i},projectedDimensions);
        Zproj = reduce(Zproj,'girard',3);
        plotFilled(Zproj,[1,2],[.8 .8 .8],'EdgeColor','none');
    end
    
    %plot initial set
    plot(options.R0,projectedDimensions,'w-','lineWidth',2);
    
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end

    %label plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end
%--------------------------------------------------------------------------
disp('comment: the rachable set becomes so thin that it is fully covered by black simulation runs.')

%example completed
completed = 1;

%------------- END OF CODE --------------
