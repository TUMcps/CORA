function completed = example_linearParam_reach_01_rlc_const()
% example_linearParam_reach_01_rlc_const - example of 
% linear parametric reachability analysis; 
% this example is also a unit test function.
%
% This example is taken from
%
% M. Althoff, B. H. Krogh, and O. Stursberg. Modeling, Design, and 
% Simulation of Systems with Uncertainties, chapter Analyzing Reachability 
% of Linear Dynamic Systems with Parametric Uncertainties, pages 69-94. 
% Springer, 2011.
%
% Syntax:  
%    example_linearParam_reach_01_rlc_const
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


%init: get matrix zonotopes of the model
[matZ_A,matZ_B] = initRLC_uTest();
matI_A = intervalMatrix(matZ_A);

%get dimension
dim=matZ_A.dim;

%compute initial set
%specify range of voltages
u0 = intervalMatrix(0,0.2);

%compute inverse of A
intA = intervalMatrix(matZ_A);
invAmid = inv(mid(intA.int)); 

%compute initial set
intB = intervalMatrix(matZ_B);
R0 = invAmid*intB*u0 + intervalMatrix(0,1e-3*ones(dim,1));

%convert initial set to zonotope
R0 = zonotope(interval(R0));

%initial set
options.x0=center(R0); %initial state for simulation
options.R0=R0; %initial state for reachability analysis

%inputs
u=intervalMatrix(1,0.01);
U = zonotope(interval(intB*u));
options.uTrans=center(U);
options.U=U+(-options.uTrans); %input for reachability analysis

%other
options.tStart=0; %start time
options.tFinal=0.7; %final time
options.intermediateOrder = 2;
options.reductionTechnique = 'girard';
options.originContained = 0;
options.timeStep = 0.002;
options.eAt = expm(matZ_A.center*options.timeStep);

options.zonotopeOrder=400; %zonotope order
options.polytopeOrder=3; %polytope order
options.taylorTerms=6;

%time step
r = options.timeStep;
maxOrder=options.taylorTerms;

%instantiate linear dynamics with constant parameters
linSys  = linParamSys(matZ_A, eye(dim), r, maxOrder);
linSys2 = linParamSys(matI_A, eye(dim), r, maxOrder);

%reachable set computations
tic
Rcont = reach(linSys, options);
tComp = toc;
disp(['computation time of reachable set using matrix zonotopes: ',num2str(tComp)]);

tic
Rcont2 = reach(linSys2, options);
tComp = toc;
disp(['computation time of reachable set using interval matrices: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(linSys2, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plot results--------------------------------------------------------------
for plotRun=1:2
    % plot different projections
    if plotRun==1
        projectedDimensions=[1 21];
    else
        projectedDimensions=[20 40];
    end
    
    figure;
    hold on
    
    %plot reachable sets
    for i=1:length(Rcont2)
        Zproj = project(Rcont2{i},projectedDimensions);
        Zproj = reduce(Zproj,'girard',3);
        plotFilled(Zproj,[1 2],[.675 .675 .675],'EdgeColor','none');
    end

    for i=1:length(Rcont)
        Zproj = project(Rcont{i},projectedDimensions);
        Zproj = reduce(Zproj,'girard',3);
        plotFilled(Zproj,[1 2],[.75 .75 .75],'EdgeColor','none');
    end

    %plot initial set
    plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
       
    %plot simulation results      
    for i=1:length(simRes.t)
        plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
    end

    %öabel plot
    xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
    ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
end

%plot results over time

figure;
hold on

%plot time elapse
for i=1:length(Rcont2)
    %get Uout 
    Uout1 = interval(project(Rcont{i},0.5*dim));
    Uout2 = interval(project(Rcont2{i},0.5*dim));
    %obtain times
    t1 = (i-1)*options.timeStep;
    t2 = i*options.timeStep;
    %generate plot areas as interval hulls
    IH1 = interval([t1; infimum(Uout1)], [t2; supremum(Uout1)]);
    IH2 = interval([t1; infimum(Uout2)], [t2; supremum(Uout2)]);

    plotFilled(IH2,[1 2],[.675 .675 .675],'EdgeColor','none');
    plotFilled(IH1,[1 2],[.75 .75 .75],'EdgeColor','none');
end

%plot simulation results
for i=1:(length(simRes.t))
    plot(simRes.t{i},simRes.x{i}(:,0.5*dim),'Color',0*[1 1 1]);
end
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------