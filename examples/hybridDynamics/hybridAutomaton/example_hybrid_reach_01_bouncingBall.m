function res = example_hybrid_reach_01_bouncingBall
% example_hybrid_reach_01_bouncingBall - example for hybrid dynamics
%    Checks the solution of the hybrid system class for the classical
%    bouncing ball example.
%
% Syntax:  
%    example_hybrid_reach_01_bouncingBall
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean  

% Author:       Matthias Althoff
% Written:      27-July-2016
% Last update:  23-December-2019
% Last revision:---

%------------- BEGIN CODE --------------


% Parameter ---------------------------------------------------------------

% problem description
params.R0 = zonotope([1;0],diag([0.05,0.05]));      % initial set
params.startLoc = 1;                                % initial location
params.tFinal = 1.7;                                % final time


% Reachability Options ----------------------------------------------------

% settings for continuous reachability 
options.timeStep = 0.05;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% settings for hybrid systems
options.guardIntersect = 'polytope';
options.enclose = {'box'}; 


% Hybrid Automaton --------------------------------------------------------

% continuous dynamics 
A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];
linSys = linearSys('linearSys',A,B,c);

% system parameters
alpha = -0.75;                  % rebound factor

% invariant set 
inv = mptPolytope([-1,0],0);

% guard sets
guard = conHyperplane([1,0],0,[0,1],0);

% reset function
reset.A = [0, 0; 0, alpha]; reset.c = zeros(2,1);

% transitions
trans{1} = transition(guard,reset,1);

% location object
loc{1} = location('loc1',inv,trans,linSys); 

% hybrid automata
HA = hybridAutomaton(loc);


% Reachability Analysis ---------------------------------------------------

tic;
R = reach(HA,params,options);
tComp = toc;

disp(['Computation time for reachable set: ',num2str(tComp),' s']);


% Simulation --------------------------------------------------------------

simRes = simulateRandom(HA,params); 



% Visualization -----------------------------------------------------------

figure; hold on;

% plot reachable set
plot(R,[1,2]);

% plot initial set
plot(params.R0,[1,2],'FaceColor','w','EdgeColor','k');

% plot simulated trajectories
plot(simRes,[1,2]);

axis([0,1.2,-6,4]);


res = true;


%------------- END OF CODE --------------
