function res = example_parallel_hybrid_01_bouncing_lever()
% example_parallel_hybrid_01_bouncing_lever - example for reachability of a
% parallel hybrid automaton. The system describes a train who's
% acceleration is controlled by a bouncing ball
%
% Syntax:  
%    example_parallel_hybrid_01_bouncing_lever
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean, true if completed

% Author:       Johann Schöpfer, Niklas Kochdumper
% Written:      06-July-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% parameters of dynamics
train_maxspeed = 2;
train_minspeed = 1.8;
train_leverpower = 1;
train_brakepower = -0.6;
ball_slowdown = -0.9;
gravity = -5;
% simTime = 10;
simTime = 3;
diff = 0.001;





% Component Train ---------------------------------------------------------

% x = [train_pos ; train_speed]
% u = [lever_pos]

% State 1: train unbraking
A = [0,1;0,0];
B = [0;train_leverpower];
c = [0;0];
C = [1,0];
flow1 = linearSys('linearSys', A, B, c, C);

inv1 = interval([-inf;-inf],[inf;train_maxspeed]);

% State 2: train braking
A = [0,1;0,train_brakepower];
B = [0;train_leverpower];
c = [0;0];
C = [1,0];
flow2 = linearSys('linearSys', A, B, c, C);

inv2 = interval([-inf;train_minspeed],[inf;inf]);

% Transitions
% reset does nothing
reset.A = eye(2);
reset.b = [0;0];

% engage brake
guard1 = halfspace([0, -1],-train_maxspeed);
trans1{1} = transition(guard1,reset,2,'','');

%release brake
guard2 = halfspace([0, 1],train_minspeed);
trans2{1} = transition(guard2,reset,1,'','');

loc_train{1} = location('nobrake',1,inv1,trans1,flow1);
loc_train{2} = location('brake',2,inv2,trans2,flow2);

HA_train = hybridAutomaton(loc_train);





% Component Bouncing Lever ------------------------------------------------

% x = [lever_pos ; lever_speed]
% u = [gravity]

% State 1: always
A = [0,1;0,0];
B = [0;1];
c = [0;0];
C = [1,0];
flow = linearSys('linearSys', A, B, c, C);

inv = interval([-diff;-inf],[inf;inf]);

% Transitions
% bounce
reset.A = [1,0;0,ball_slowdown];
reset.b = [0;0];
temp = halfspace([1,0],0); 
guard = constrainedHyperplane(temp,[0 1],0);
trans{1} = transition(guard,reset,1,'','');

loc_ball{1} = location('always',1,inv,trans,flow);

HA_ball = hybridAutomaton(loc_ball);




% Composition -------------------------------------------------------------

comp{1} = HA_train;
comp{2} = HA_ball;

% Composition of states
stateBinds{1} = [1;2];
stateBinds{2} = [3;4];

% Connection of the subcomponents
inputBinds{1} = [2,1];
inputBinds{2} = [0,1];

PHA = parallelHybridAutomaton(comp,stateBinds,inputBinds);




% Simulation --------------------------------------------------------------

% general options
options.tStart = 0; %start time
options.tFinal = simTime; %final time
options.taylorTerms = 10;
options.zonotopeOrder = 30;
options.polytopeOrder = 10;
options.reductionTechnique = 'girard';
options.isHyperplaneMap = 0;
options.guardIntersect = 'zonoGirard';
options.enclosureEnables = [1,2];
options.originContained = 0;
options.linAlg = 2;

% PHA options
options_pha = options;
options_pha.x0 = [0; 0; 2; 0]; %initial state for simulation
options_pha.R0 = zonotope([options_pha.x0, diag([0.01 0.01 0.01 0.01])]);
options_pha.startLoc = {1,1}; %initial location
options_pha.finalLoc = {0,0}; %0: no final location

% PHA input
% Explanation of component-wise input:
% Input domain of PHA: U = [u1 u2 u3 u4 u5]
% inputCompMap is a mapping: [1,numInputs] -> [0,numComponents]
% (example: inputCompMap(1:5) = [1  1  3  2  3 ], for a 3-component PHA)
% => Input domain of Component c: U_c = [u_i | m(i) == c]
% => Additional global input domain: U_glob = [u_i | m(i) == 0]
% (example: U_1 = [u1 u2]; U_2 = [u4]; U_3 = [u3 u5]; U_glob = [])
% => for all component-locations l: uCompLoc{c}{l} must specify U_c
% => uGlob must specify U_glob
% (example: uCompLoc{1}{*} ∈ R^2; uCompLoc{3}{*} ∈ R^2; uCompLoc{2}{*} ∈ R
%           uGlob ∈ R^0)

options_pha.inputCompMap = [2];
options_pha.uCompLoc{2}{1} = gravity;
options_pha.UCompLoc{2}{1} = zonotope(0);
options_pha.uCompLocTrans{2}{1} = gravity;
options_pha.uGlobTrans = [];
options_pha.uGlob = [];
options_pha.UGlob = [];

PHA = simulateRandom(PHA,10,0.5,0.5,20,options_pha);




% Reachability Analysis ---------------------------------------------------

options_pha.UCompLoc{2}{1} = zonotope(0);
options_pha.uCompLocTrans{2}{1} = gravity;
options_pha.uGlobTrans = [];

options_pha.timeStep = 0.01;

PHA = reach(PHA, options_pha);




% Simulation Bouncing Lever -----------------------------------------------

options.uLoc{1} = gravity;
options.uLocTrans{1} = options.uLoc{1}; %input center for reachability analysis
options.Uloc{1} = zonotope(0); %input deviation for reachability analysis
options.finalLoc = -1;
options.startLoc = 1;
options.x0 = options_pha.x0(3:4);
options.projectedDimensions = [1,2];
options.plotType = 'g';

%simulate hybrid automaton
HA_ball = simulate(HA_ball,options); 



% Visualization -----------------------------------------------------------

% plot train model
figure 
box on
options_pha.projectedDimensions = [1 2];
options_pha.plotType = 'r';
plot(PHA,'reachableSet',options_pha);
options_pha.plotType = 'k';
plot(PHA,'simulation',options_pha);
xlabel('$x_{train}$','interpreter','latex','FontSize',20);
ylabel('$v_{train}$','interpreter','latex','FontSize',20);
title('Train');

% plot the bouncing ball model
figure 
box on
options_pha.projectedDimensions = [3 4];
% options_pha.plotType = 'xr';
options_pha.plotType = 'g';
plot(PHA,'reachableSet',options_pha);
options_pha.plotType = 'k';
plot(PHA,'simulation',options_pha); %plot simulation
% plot(HA_ball,'simulation',options);
xlabel('$x_{ball}$','interpreter','latex','FontSize',20);
ylabel('$v_{ball}$','interpreter','latex','FontSize',20);
title('Bouncing Ball');

% plot the combination
figure 
box on
options_pha.projectedDimensions = [1 3];
options_pha.plotType = 'b';
plot(PHA,'reachableSet',options_pha);
options_pha.plotType = 'k';
plot(PHA,'simulation',options_pha); %plot simulation
xlabel('$x_{train}$','interpreter','latex','FontSize',20);
ylabel('$x_{ball}$','interpreter','latex','FontSize',20);

res = 1;

end

%------------- END OF CODE --------------