function completed = example_nonlinear_reach_03_vanDerPol()
% example_nonlinear_reach_03_vanDerPol - example of nonlinear reachability 
% analysis; this example is also a unit test function.
%
% This example can be found in Sec. 3.4.5 of
% M. Althoff, “Reachability analysis and its application to the safety 
% assessment of autonomous cars�?, Dissertation, Technische Universität 
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
% A new technique for computing this example with less spliiting has bee
% published in
%
% M. Althoff. Reachability analysis of nonlinear systems using conservative 
% polynomialization and non-convex sets. In Hybrid Systems: Computation 
% and Control, pages 173-182, 2013
%
% Syntax:  
%    example_nonlinear_reach_03_vanDerPol()
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
% Written:      26-June-2009
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


%set options---------------------------------------------------------------
options.tStart=0; %start time
%options.tFinal=0.5; %final time
options.tFinal=6.74; %final time
options.x0=[1.4; 2.3]; %initial state for simulation

Z0{1}=zonotope([1.4 0.3 0; 2.3 0 0.05]); %initial state for reachability analysis
%Z0{1}=zonotope([1.4 0.6 0; 2.3 0 0.05]); %initial state for reachability analysis
options.R0 = zonotopeBundle(Z0);

options.timeStep=0.02; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.intermediateOrder = 10;
options.errorOrder = 5;
options.polytopeOrder=1.5; %polytope order

options.advancedLinErrorComp = 1;
options.tensorOrder=2;
options.reductionTechnique='girard';
options.filterLength = [10, 5];
options.maxError=0.05*[1; 1];
options.reductionInterval=100;
options.verbose = 1;

%uncertain inputs
options.uTrans = 0;
options.U = zonotope(0); %input for reachability analysis
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
vanderPol=nonlinearSys(2,1,@vanderPolEq,options); %initialize van-der-Pol oscillator
%--------------------------------------------------------------------------
      
tic
%compute reachable set
Rcont = reach(vanderPol, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

%create random simulations; RRTs would provide better results, but are
%computationally more demanding
runs = 60;
fractionVertices = 0.5;
fractionInputVertices = 0.5;
inputChanges = 6;
simRes = simulate_random(vanderPol, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plotting message
disp('Start plotting; takes a while since plotting acceleration for zonotope bundles not yet implemented');

%plot results--------------------------------------------------------------
projectedDimensions=[1 2];
plotOrder = 20;
    
figure;
hold on

%plot reachable sets 
for i=1:length(Rcont)
    for j=1:length(Rcont{i})
        Zproj = reduce(Rcont{i}{j},'girard',plotOrder);
        plotFilled(Zproj,projectedDimensions,[.8 .8 .8],'EdgeColor','none');
    end
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
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------