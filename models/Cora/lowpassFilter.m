function pHA = lowpassFilter()
% lowpassFilter - creates an instance of a parallel hybrid automaton from
%    the model of a lowpass filter with 2 components and 3 locations each
%
% Syntax:  
%    PHA = lowpassFilter()
%
% Inputs:
%    -
%
% Outputs:
%    pHA - object of parallelHybridAutomaton class

% Author:       Mark Wetzlinger
% Written:      16-December-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Low-Pass-Filter 1 -------------------------------------------------------

% Location 1
A11 = [  -1.021587e+01  0.000000e+00 ; 
             0.000000e+00 -2.002418e+02] ;  
B11 = [0;1] ; 
c11 = [20.13684457; -400.4898478];
C11 = [1,0];

linSys11  =  linearSys('linearSys',A11 ,B11, c11, C11); 

inv = interval(  [  -1.097840e+00 ; -2.698928e+00 ; ] ,  [   1.097840e+00 ;  2.698928e+00 ; ]  ) ; 

guard1   = interval(  [  -1.107840e+00 ; -2.698928e+00 ; ] ,  [  -1.087840e+00 ;  2.698928e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.c =   [   1.071580e+00 ; -1.301798e+00 ; ] ; 
tran = transition( guard1, reset,2);

guard2   = interval(  [   1.087840e+00 ; -2.698928e+00 ; ] ,  [   1.107840e+00 ;  2.698928e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.c =   [  -1.071580e+00 ;  1.301798e+00 ; ] ; 
tran(2) = transition( guard2, reset,3) ;

loc(1)  =  location('loc11',inv,tran,linSys11); 

% Location 2
A21 = [  -1.458328e+01  0.000000e+00 ; 
             0.000000e+00 -2.053664e+02] ;  
B21 = [0;1] ;  
c21 = [33.95001456; -632.6889324];
C21 = [0,1];

linSys21  =  linearSys('linearSys',A21 ,B21, c21, C21); 

inv = interval(  [  -1.629461e+00 ; -2.653899e+00 ; ] ,  [   3.718705e-01 ;  2.705073e+00 ; ]  ) ; 
 
guard1   = interval(  [   3.618705e-01 ; -2.653899e+00 ; ] ,  [   3.818705e-01 ;  2.705073e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.c =   [  -0.000000e+00 ; -0.000000e+00 ; ] ; 
tran = transition( guard1, reset,1);

loc(2)  =  location('loc21',inv,tran,linSys21);

% Location 3
A31 = [  -1.458328e+01  0.000000e+00 ; 
             0.000000e+00 -2.053664e+02] ;  
B31 = [0;1] ;  
c31 = [7.87473456; -146.7527324];
C31 = [0,1];

linSys31  =  linearSys('linearSys',A31 ,B31, c31, C31); 

inv = interval(  [  -3.718705e-01 ; -2.705073e+00 ; ] ,  [   1.629461e+00 ;  2.653899e+00 ; ]  ) ; 

guard1   = interval(  [  -3.818705e-01 ; -2.705073e+00 ; ] ,  [  -3.618705e-01 ;  2.653899e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.c =   [  -0.000000e+00 ; -0.000000e+00 ; ] ; 
tran = transition( guard1, reset,1);

loc(3)  =  location('loc31',inv,tran,linSys31); 

% Hybrid automaton
HA_1 = hybridAutomaton(loc); 



% Low-Pass-Filter 2 -------------------------------------------------------

% Location 1
A11 = [  -1.265885e+01  0.000000e+00 ; 
             0.000000e+00 -4.725145e+01] ;  
B11 = [0;1] ;
c11 = [-2.92283193; 8.0628348];
C11 = [0,1];

linSys11  =  linearSys('linearSys',A11 ,B11, c11, C11);

inv = interval(  [  -1.350712e+00 ; -2.504938e+00 ; ] ,  [   9.670362e-01 ;  1.562432e+00 ; ]  ) ; 

guard1   = interval(  [   9.570362e-01 ; -2.504938e+00 ; ] ,  [   9.770362e-01 ;  1.562432e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.c =   [  -0.000000e+00 ; -0.000000e+00 ; ] ; 
tran = transition( guard1, reset,3);

loc(1)  =  location('loc11',inv,tran,linSys11); 

% Location 2
A21 = [  -1.265885e+01  0.000000e+00 ; 
             0.000000e+00 -4.725145e+01] ;  
B21 = [0;1] ;  
c21 = [-43.04013293; 118.7291948];
C21 = [0,1];

linSys21  =  linearSys('linearSys',A21 ,B21, c21, C21); 

inv = interval(  [  -9.670362e-01 ; -1.562432e+00 ; ] ,  [   1.350712e+00 ;  2.504938e+00 ; ]  ) ; 
 
guard1   = interval(  [  -9.770362e-01 ; -1.562432e+00 ; ] ,  [  -9.570362e-01 ;  2.504938e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.c =   [  -0.000000e+00 ; -0.000000e+00 ; ] ; 
tran = transition( guard1, reset,3);

loc(2)  =  location('loc21',inv,tran,linSys21);  

% Location 3
A31 = [  -1.007836e+01  0.000000e+00 ; 
             0.000000e+00 -4.016005e+01] ;  
B31 = [0;1] ;  
c31 = [-26.0284288; 79.7753009];
C31 = [0,1];

linSys31  =  linearSys('linearSys',A31 ,B31, c31, C31); 

inv = interval(  [  -2.251641e+00 ; -2.580417e+00 ; ] ,  [   2.251641e+00 ;  2.580417e+00 ; ]  ) ; 
 
guard1   = interval(  [  -2.261641e+00 ; -2.580417e+00 ; ] ,  [  -2.241641e+00 ;  2.580417e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.c =   [   2.089598e+00 ; -1.797423e+00 ; ] ; 
tran = transition( guard1, reset,1);

guard2   = interval(  [   2.241641e+00 ; -2.580417e+00 ; ] ,  [   2.261641e+00 ;  2.580417e+00 ; ]  ) ; 
reset.A = eye(2,2); reset.c =   [  -2.089598e+00 ;  1.797423e+00 ; ] ; 
tran(2) = transition( guard2, reset,2) ;

loc(3)  =  location('loc31',inv,tran,linSys31); 
 
% Hybrid automata 
HA_2 = hybridAutomaton(loc); 





% Parallel Hybrid Automaton -----------------------------------------------

% list of subcomponents
comp(1) = HA_1;
comp(2) = HA_2;

% connection of the subcomponents 
inputBinds{1} = [0,1];
inputBinds{2} = [1,1];

% construct parallel hybrd automaton
pHA = parallelHybridAutomaton(comp,inputBinds);

end

%------------- END OF CODE --------------