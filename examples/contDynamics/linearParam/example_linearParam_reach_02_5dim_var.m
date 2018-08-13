function completed = example_linearParam_reach_02_5dim_var()
% example_linearParam_reach_02_5dim_var - example of 
% linear parametric reachability analysis where the parameters can vary 
% over time; this example is also a unit test function.
%
% This example is taken from
%
% Althoff, M.; Le Guernic, C. & Krogh, B. H. Reachable Set Computation for 
% Uncertain Time-Varying Linear Systems Hybrid Systems: Computation and 
% Control, 2011, 93-102
%
% Syntax:  
%    example_linearParam_reach_02_5dim_var
%
% Inputs:
%    none
%
% Outputs:
%    none
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      03-October-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


dim = 5;

%set options --------------------------------------------------------------
options.tStart = 0; %start time
options.tFinal = 5; %final time
options.x0 = ones(dim,1); %initial state for simulation
options.R0 = zonotope([options.x0,0.1*eye(length(options.x0))]); %initial state for reachability analysis

options.timeStep = 0.05; %time step size for reachable set computation
options.taylorTerms = 4; %number of taylor terms for reachable sets
options.intermediateOrder = 2;
options.zonotopeOrder = 20; %zonotope order
options.originContained = 0;
options.reductionTechnique = 'girard';

uTrans = zeros(dim,1);
options.uTrans = uTrans; %center of input set
options.U = zonotope([zeros(dim,1),0.1*eye(dim)]); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
Acenter = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
Arad{1} = [0.1 0.1 0 0 0; 0.1 0.1 0 0 0; 0 0 0.1 0.1 0; 0 0 0.1 0.1 0; 0 0 0 0 0.1];
matZ_A = matZonotope(Acenter,Arad);
matI_A = intervalMatrix(matZ_A);

fiveDimSys_zono = linParamSys(matZ_A, 1, options.timeStep, options.taylorTerms,'varParam'); %instantiate system
fiveDimSys_int = linParamSys(matI_A, 1, options.timeStep, options.taylorTerms,'varParam'); %instantiate system
%--------------------------------------------------------------------------

%compute reachable set using zonotopes
tic
Rcont_zono = reach(fiveDimSys_zono, options);
tComp = toc;
disp(['computation time of reachable set for zonotopic matrix: ',num2str(tComp)]);

%compute reachable set using zonotopes
tic
Rcont_int = reach(fiveDimSys_int, options);
tComp = toc;
disp(['computation time of reachable set for interval matrix: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(fiveDimSys_zono, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plot results--------------------------------------------------------------
for plotRun = 1:2
    % plot different projections
    if plotRun == 1
        projectedDimensions=[2 3];
    elseif plotRun == 2
        projectedDimensions=[4 5];     
    end 
    
    figure;
    hold on
    
    %plot reachable sets of interval matrix
    for i=1:length(Rcont_int)
        plotFilled(Rcont_int{i},projectedDimensions,[.675 .675 .675],'EdgeColor','none');
    end

    %plot reachable sets 
    for i=1:length(Rcont_zono)
        plotFilled(Rcont_zono{i},projectedDimensions,[.8 .8 .8],'EdgeColor','none');
    end
    
    %plot initial set
    plot(options.R0,projectedDimensions,'w-','lineWidth',2);
    
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