function completed = example_linear_reach_01_5dim()
% example_linear_reach_01_5dim - example of linear reachability 
% analysis with uncertain inputs; this example is also a unit test
% function.
%
% This example can be found in Sec. 3.2.3 of
% M. Althoff, “Reachability analysis and its application to the safety 
% assessment of autonomous cars”, Dissertation, Technische Universität 
% München, 2010, 
% http://nbn-resolving.de/urn/resolver.pl?urn:nbn:de:bvb:91-diss-20100715-963752-1-4
%
% Syntax:  
%    example_linear_reach_01_5dim
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
% Written:      17-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim=5;

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=5; %final time
options.x0=ones(dim,1); %initial state for simulation
options.R0=zonotope([options.x0,0.1*eye(length(options.x0))]); %initial state for reachability analysis

options.timeStep=0.04; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=200; %zonotope order
options.originContained=0;
options.reductionTechnique='girard';

uTrans=[1; 0; 0; 0.5; -0.5];
options.uTrans=uTrans; %center of input set
options.U=0.5*zonotope([zeros(5,1),diag([0.2, 0.5, 0.2, 0.5, 0.5])]); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
A=[-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B=1;
fiveDimSys=linearSys('fiveDimSys',A,B); %initialize system
%--------------------------------------------------------------------------

%compute reachable set using zonotopes
tic
Rcont = reach(fiveDimSys, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(fiveDimSys, options, runs, fractionVertices, fractionInputVertices, inputChanges);

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
        plotFilled(Rcont{i},projectedDimensions,[.8 .8 .8],'EdgeColor','none');
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
