function res = example_parallel_hybrid_02_lowPassFilter()
% example_parallel_hybrid_02_lowPassFilter - example for reachability of a
% parallel hybrid automaton. The system consists of two piecewise linear
% low-pass filters that are connected in series
%
% Syntax:  
%    example_parallel_hybrid_02_lowPassFilter
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean, true if completed

% Author:       Niklas Kochdumper
% Written:      06-July-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------




% Low-Pass-Filter 1 -------------------------------------------------------

% Invariant 1
A11 = [  -1.021587e+01  0.000000e+00 ; 
             0.000000e+00 -2.002418e+02] ;  
B11 = [0;1] ; 
c11 = [20.13684457; -400.4898478];
C11 = [1,0];
D11 = 0;
k11 = 0;

linSys11  =  linearSys('linearSys',A11 ,B11, c11, C11, D11, k11); 

inv = interval(  [  -1.097840e+00 ; -2.698928e+00 ; ] ,  [   1.097840e+00 ;  2.698928e+00 ; ]  ) ; 

guard1   = interval(  [  -1.107840e+00 ; -2.698928e+00 ; ] ,  [  -1.087840e+00 ;  2.698928e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.b =   [   1.071580e+00 ; -1.301798e+00 ; ] ; 
tran = {transition( guard1, reset,21,'u21','y21')};

guard2   = interval(  [   1.087840e+00 ; -2.698928e+00 ; ] ,  [   1.107840e+00 ;  2.698928e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.b =   [  -1.071580e+00 ;  1.301798e+00 ; ] ; 
tran{ 2 } = transition( guard2, reset,31,'u31','y31') ;

loc{11}  =  location('loc11',11,inv,tran,linSys11); 

% figure
% hold on
% plot(inv,[1 2],'color',[5.878043e-01 8.290435e-01 1.777418e-01]);  
% plot(guard1,[1 2],'color',[2.513277e-01 1.484674e-01 5.216626e-01]); 
% plot(guard2,[1 2],'color',[2.633118e-01 8.782764e-01 9.981223e-01]); 



% Invariant 2
A21 = [  -1.458328e+01  0.000000e+00 ; 
             0.000000e+00 -2.053664e+02] ;  
B21 = [0;1] ;  
c21 = [33.95001456; -632.6889324];
C21 = [0,1];
D21 = 0;
k21 = 0;

linSys21  =  linearSys('linearSys',A21 ,B21, c21, C21, D21, k21); 

inv = interval(  [  -1.629461e+00 ; -2.653899e+00 ; ] ,  [   3.718705e-01 ;  2.705073e+00 ; ]  ) ; 
 
guard1   = interval(  [   3.618705e-01 ; -2.653899e+00 ; ] ,  [   3.818705e-01 ;  2.705073e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.b =   [  -0.000000e+00 ; -0.000000e+00 ; ] ; 
tran = {transition( guard1, reset,11,'u11','y11')};

loc{21}  =  location('loc21',21,inv,tran,linSys21); 

% plot(inv,[1 2],'color',[8.280972e-01 9.166018e-01 2.338734e-01]);
% plot(guard1,[1 2],'color',[8.480908e-01 7.699192e-01 8.049495e-01]); 



% Invariant 3
A31 = [  -1.458328e+01  0.000000e+00 ; 
             0.000000e+00 -2.053664e+02] ;  
B31 = [0;1] ;  
c31 = [7.87473456; -146.7527324];
C31 = [0,1];
D31 = 0;
k31 = 0;

linSys31  =  linearSys('linearSys',A31 ,B31, c31, C31, D31, k31); 

inv = interval(  [  -3.718705e-01 ; -2.705073e+00 ; ] ,  [   1.629461e+00 ;  2.653899e+00 ; ]  ) ; 

guard1   = interval(  [  -3.818705e-01 ; -2.705073e+00 ; ] ,  [  -3.618705e-01 ;  2.653899e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.b =   [  -0.000000e+00 ; -0.000000e+00 ; ] ; 
tran = {transition( guard1, reset,11,'u11','y11')};

loc{31}  =  location('loc31',31,inv,tran,linSys31); 

% plot(inv,[1 2],'color',[5.287919e-02 6.497299e-01 2.751759e-02]);  
% plot(guard1,[1 2],'color',[2.824982e-01 1.550679e-01 3.950005e-01]); 

% Hybrid automata 
HA_1 = hybridAutomaton(loc); 






% Low-Pass-Filter 2 -------------------------------------------------------

% Invariant 1
A11 = [  -1.265885e+01  0.000000e+00 ; 
             0.000000e+00 -4.725145e+01] ;  
B11 = [0;1] ;
c11 = [-2.92283193; 8.0628348];
C11 = [0,1];
D11 = 0;
k11 = 0;

linSys11  =  linearSys('linearSys',A11 ,B11, c11, C11, D11, k11);

inv = interval(  [  -1.350712e+00 ; -2.504938e+00 ; ] ,  [   9.670362e-01 ;  1.562432e+00 ; ]  ) ; 

guard1   = interval(  [   9.570362e-01 ; -2.504938e+00 ; ] ,  [   9.770362e-01 ;  1.562432e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.b =   [  -0.000000e+00 ; -0.000000e+00 ; ] ; 
tran = {transition( guard1, reset,31,'u31','y31')};

loc{11}  =  location('loc11',11,inv,tran,linSys11); 

% figure
% hold on
% plot(inv,[1 2],'color',[8.290435e-01 1.777418e-01 2.513277e-01]);  
% plot(guard1,[1 2],'color',[1.484674e-01 5.216626e-01 2.633118e-01]); 



% Invariant 2
A21 = [  -1.265885e+01  0.000000e+00 ; 
             0.000000e+00 -4.725145e+01] ;  
B21 = [0;1] ;  
c21 = [-43.04013293; 118.7291948];
C21 = [0,1];
D21 = 0;
k21 = 0;

linSys21  =  linearSys('linearSys',A21 ,B21, c21, C21, D21, k21); 

inv = interval(  [  -9.670362e-01 ; -1.562432e+00 ; ] ,  [   1.350712e+00 ;  2.504938e+00 ; ]  ) ; 
 
guard1   = interval(  [  -9.770362e-01 ; -1.562432e+00 ; ] ,  [  -9.570362e-01 ;  2.504938e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.b =   [  -0.000000e+00 ; -0.000000e+00 ; ] ; 
tran = {transition( guard1, reset,31,'u31','y31')};

loc{21}  =  location('loc21',21,inv,tran,linSys21); 

% plot(inv,[1 2],'color',[8.782764e-01 9.981223e-01 8.280972e-01]); 
% plot(guard1,[1 2],'color',[9.166018e-01 2.338734e-01 8.480908e-01]); 


% Invariant 3
A31 = [  -1.007836e+01  0.000000e+00 ; 
             0.000000e+00 -4.016005e+01] ;  
B31 = [0;1] ;  
c31 = [-26.0284288; 79.7753009];
C31 = [0,1];
D31 = 0;
k31 = 0;

linSys31  =  linearSys('linearSys',A31 ,B31, c31, C31, D31, k31); 

inv = interval(  [  -2.251641e+00 ; -2.580417e+00 ; ] ,  [   2.251641e+00 ;  2.580417e+00 ; ]  ) ; 
 
guard1   = interval(  [  -2.261641e+00 ; -2.580417e+00 ; ] ,  [  -2.241641e+00 ;  2.580417e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.b =   [   2.089598e+00 ; -1.797423e+00 ; ] ; 
tran = {transition( guard1, reset,11,'u11','y11')};

guard2   = interval(  [   2.241641e+00 ; -2.580417e+00 ; ] ,  [   2.261641e+00 ;  2.580417e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.b =   [  -2.089598e+00 ;  1.797423e+00 ; ] ; 
tran{ 2 } = transition( guard2, reset,21,'u21','y21') ;

loc{31}  =  location('loc31',31,inv,tran,linSys31); 

% plot(inv,[1 2],'color',[7.699192e-01 8.049495e-01 5.287919e-02]); 
% plot(guard1,[1 2],'color',[6.497299e-01 2.751759e-02 2.824982e-01]); 
% plot(guard2,[1 2],'color',[1.550679e-01 3.950005e-01 7.855495e-01]); 

 
%Hybrid automata 
HA_2 = hybridAutomaton(loc); 





% Parallel Hybrid Automaton -----------------------------------------------

comp{1} = HA_1;
comp{2} = HA_2;

% composition of the overall state vector
stateBinds{1} = [1;2];
stateBinds{2} = [3;4];

% connection of the subcomponents 
inputBinds{1} = [0,1];
inputBinds{2} = [1,1];

PHA = parallelHybridAutomaton(comp,stateBinds,inputBinds);




% Options -----------------------------------------------------------------

% General Options
options.tStart = 0; 
% options.tFinal = 0.8;
options.tFinal = 0.4;
options.taylorTerms = 8; 
options.zonotopeOrder = 9; 
options.polytopeOrder = 10; 
options.errorOrder = 2; 
options.originContained = 1; 
options.reductionTechnique = 'girard'; 
options.enclosureEnables = [1]; 
options.isHyperplaneMap = 0; 
% options.guardIntersect = 'polytope';
options.guardIntersect = 'polytope';

% Options for hybrid automata
options.x0 = [0;0;0;0] ; 
options.R0 = zonotope([options.x0,diag([0.01,0.01,0.1,0.1])]);  
options.timeStep = 1e-04 ; 
options.startLoc = {11, 31}; 
options.finalLoc = {0, 0}; 

% System Inputs
options.inputCompMap = [0];
options.uCompLoc = [];
options.UCompLoc = [];
options.uCompLocTrans = [];
options.uGlobTrans = 0;
options.uGlob = 0;
options.UGlob = zonotope(0);




% Simulation --------------------------------------------------------------

PHA = simulateRandom(PHA,10,0.5,0.5,20,options);




% Reachability Analysis ---------------------------------------------------

PHA = reach(PHA,options);




% Visualization -----------------------------------------------------------

% plot filter 1
figure 
box on
hold on
options.projectedDimensions = [1 2];
options.plotType = 'r';
plot(PHA,'reachableSet',options);
options.plotType = 'k';
plot(PHA,'simulation',options);
xlabel('$x_{1}$','interpreter','latex','FontSize',20);
ylabel('$x_{2}$','interpreter','latex','FontSize',20);
title('Filter 1');

% plot filter 2
figure 
box on
options.projectedDimensions = [3 4];
options.plotType = 'g';
plot(PHA,'reachableSet',options);
options.plotType = 'k';
plot(PHA,'simulation',options);
xlabel('$x_{1}$','interpreter','latex','FontSize',20);
ylabel('$x_{2}$','interpreter','latex','FontSize',20);
title('Filter 2');



% Visited Locations -------------------------------------------------------

visitedLocations = get(PHA,'locationReach');

res = 1;

%------------- END OF CODE --------------