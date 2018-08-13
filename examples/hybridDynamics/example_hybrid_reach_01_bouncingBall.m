function completed = example_hybrid_reach_01_bouncingBall
% example_hybrid_reach_01_bouncingBall - example for hybrid
% dynamics; this example is also a unit test function.
%
% Checks the solution of the hybrid system class for the classical bouncing 
% ball example.
%
% Syntax:  
%    example_hybrid_reach_01_bouncingBall
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
% Written:      27-July-2016
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
%specify linear system of bouncing ball
A = [0 1; 0 0];
B = [0; 0]; % no inputs
c = [0; -9.81]; % constant part of affine dynamics
linSys = linearSys('linearSys',A,B,c);

%define large and small distance
dist = 1e3;
eps = 1e-6;
alpha = -0.75; %rebound factor

%invariant
inv = interval([-2*eps; -dist], [dist; dist]);
%guard sets
guard = interval([-eps; -dist], [0; -eps]); 
%resets
reset.A = [0, 0; 0, alpha]; reset.b = zeros(2,1);
%transitions
trans{1} = transition(guard,reset,1,'a','b'); %--> next loc: 1; 'a', 'b' are dummies
%specify location
loc{1} = location('loc1',1,inv,trans,linSys); 
%specify hybrid automata
HA = hybridAutomaton(loc); % for "geometric intersection"
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
