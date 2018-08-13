function completed = example_hybrid_reach_01_bouncingBall_converted
% example_hybrid_reach_01_bouncingBall_converted - example for hybrid
% dynamics; this example shows how to include automatically generated
% models
%
% Syntax:  
%    example_hybrid_reach_01_bouncingBall_converted
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
% Written:      07-August-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


%set options---------------------------------------------------------------
options.x0 = [1; 0]; %initial state for simulation
options.R0 = zonotope([options.x0, diag([0.05, 0.05])]); %initial state for reachability analysis
options.startLoc = 1; %initial location
options.finalLoc = 0; %0: no final location
options.tStart = 0; %start time
options.tFinal = 1.7; %final time
options.timeStepLoc{1} = 0.05; %time step size for reachable set computation in location 1
options.taylorTerms = 10;
options.zonotopeOrder = 20;
options.polytopeOrder = 10;
options.errorOrder=2;
options.reductionTechnique = 'girard';
options.isHyperplaneMap = 0;
options.guardIntersect = 'polytope';
options.enclosureEnables = 5; %choose enclosure method(s)
options.originContained = 0;
%--------------------------------------------------------------------------


%specify hybrid automaton--------------------------------------------------
% converetd hybrid automaton model of the bouncing ball obtained from 
% "spaceex2cora(’bball.xml’);"
HA = bball; 
%--------------------------------------------------------------------------

%set input:
options.uLoc{1} = 0; % no inputs
options.uLocTrans{1} = 0; % no inputs
options.Uloc{1} = zonotope(0); % no inputs

%simulate hybrid automaton
HA = simulate(HA,options); 

%compute reachable set
[HA] = reach(HA,options);

%choose projection and plot------------------------------------------------
figure 
hold on
options.projectedDimensions = [1 2];
options.plotType = 'b';
plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA,'simulation',options); %plot simulation
axis([0,1.2,-6,4]);
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------
