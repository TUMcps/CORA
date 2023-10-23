function HA = rendezvous_nonlinear_passive_nondet(~)
% rendezvous_nonlinear_passive_nondet - nonlinear spacecraft-rendezvous
%    system with nondeterministic switching (see Sec. 3.6 in [1])
%
% Syntax:  
%    HA = rendezvous_nonlinear_passive_nondet()
%
% Inputs:
%    ---
%
% Outputs:
%    HA - hybridAutomaton object
% 
% References:
%    [1] L. Geretti, “ARCH-COMP20 Category Report: Continuous and Hybrid 
%        Systems with Nonlinear Dynamics", 2020

% Author:        Niklas Kochdumper
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------


%% Generated on 04-Jun-2018

%---------------Automaton created from Component 'system'------------------

%% Interface Specification:
%   This section clarifies the meaning of state & input dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (system.SpaceCraft):
%  state x := [x; y; vx; vy; t]
%  input u := [uDummy]

%----------------------Component system.SpaceCraft-------------------------

%-------------------------------State P1-----------------------------------

%% equation:
%   
%               x'==vx &
%               y'==vy &
%               vx'== (n^2 + K1_11/m_c)*x + (2*n + K1_14/m_c)*vy + K1_12/m_c * y + K1_13/m_c * vx + mu/r^2 - mu/sqrt((r+x)^2-y^2)^3 * (r+x)&
%               vy'== (n^2 + K1_22/m_c)*y + (K1_23/m_c -2*n)*vx + K1_21/m_c * x + K1_24/m_c * vy - mu/sqrt((r+x)^2-y^2)^3 * y &
%               t'==1
%          
dynamics = nonlinearSys(@rendezvous_nonlinear_passive_St1_FlowEq); 

%% equation:
%   t<=150 & x<=-100
invA = ...
[0,0,0,0,1;1,0,0,0,0];
invb = ...
[150;-100];
invOpt = struct('A', invA, 'b', invb);
inv = polytope(invA, invb);

trans = {};
%% equation:
%   
resetA = ...
[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetc = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'c', resetc);

%% equation:
%   y>=-100 & x+y >=-141.1 & x>=-100 & y-x<=141.1 & y<=100 & x+y<=141.1 & x<=100 & y-x>=-141.1 & t < 120
guard = conHyperplane([-1,0,0,0,0],100);

trans = transition(guard, reset, 2);

%% equation:
%   
resetA = ...
[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetc = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'c', resetc);

%% equation:
%   t >= 120 & t <= 150
guard = polytope([0 0 0 0 -1; 0 0 0 0 1],[-120; 150]);

trans(2) = transition(guard, reset, 3);

loc(1) = location('S1', inv, trans, dynamics);



%-------------------------------State P2-----------------------------------

%% equation:
%   
%               x'==vx &
%               y'==vy &
%               vx'== (n^2 + K2_11/m_c)*x + (2*n + K2_14/m_c)*vy + K2_12/m_c * y + K2_13/m_c * vx + mu/r^2 - mu/sqrt((r+x)^2-y^2)^3 * (r+x)&
%               vy'== (n^2 + K2_22/m_c)*y + (K2_23/m_c -2*n)*vx + K2_21/m_c * x + K2_24/m_c * vy - mu/sqrt((r+x)^2-y^2)^3 * y &
%               t'==1
%          
dynamics = nonlinearSys(@rendezvous_nonlinear_passive_St2_FlowEq); 

%% equation:
%   t<150 & y>=-100 & x+y>=-141.1 & x>=-100 & y-x<=141.1 & y<=100 & x+y<=141.1 & x<=100 & y-x>=-141.1
invA = ...
[0,0,0,0,1;0,-1,0,0,0;-1,-1,0,0,0;-1,0,0,0,0;-1,1,0,0,0;0,1,0,0,0;1,1,0,...
0,0;1,0,0,0,0;1,-1,0,0,0];
invb = ...
[150;100;141.1;100;141.1;100;141.1;100;141.1];
invOpt = struct('A', invA, 'b', invb);
inv = polytope(invA, invb);

trans = transition();
%% equation:
%   
resetA = ...
[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetc = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'c', resetc);

%% equation:
%   t >= 120 & t <= 150
guard = polytope([0 0 0 0 -1; 0 0 0 0 1],[-120; 150]);

trans = transition(guard, reset, 3);

loc(2) = location('S2', inv, trans, dynamics);



%-----------------------------State Passive--------------------------------

%% equation:
%   
%               x'==vx & 
%               y'==vy & 
%               vx'== n^2 * x + 2*n*vy + mu/r^2 - mu/sqrt((r+x)^2-y^2)^3 * (r+x) &
%               vy'== n^2*y - 2*n*vx - mu/sqrt((r+x)^2-y^2)^3 * y &
%               t'==1
%         
dynamics = nonlinearSys(@rendezvous_nonlinear_passive_St3_FlowEq); 

%% equation:
%   x<=10000 & x>=-10000
invA = ...
[1,0,0,0,0;-1,0,0,0,0];
invb = ...
[10000;10000];
invOpt = struct('A', invA, 'b', invb);
inv = polytope(invA, invb);

trans = transition();
loc(3) = location('S3', inv, trans, dynamics);



HA = hybridAutomaton(loc);


end

%------------- END OF CODE --------------