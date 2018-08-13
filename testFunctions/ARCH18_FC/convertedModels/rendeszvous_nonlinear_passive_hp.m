function HA = rendeszvous_nonlinear_passive_hp(~)


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
dynOpt = struct('tensorOrder',1);
dynamics = nonlinearSys(5,1,@rendeszvous_nonlinear_passive_St1_FlowEq,dynOpt); 

%% equation:
%   t<=125&x<=-100
invA = ...
[0,0,0,0,1;1,0,0,0,0];
invb = ...
[125;-100];
invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

trans = {};
%% equation:
%   
resetA = ...
[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetb = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   y>=-100 & x+y >=-141.1 & x>=-100 & y-x<=141.1 & y<=100 & x+y<=141.1 & x<=100 & y-x>=-141.1 & t < 120
guard = halfspace([-1,0,0,0,0],100);

trans{1} = transition(guard, reset, 2, 'dummy', 'names');

%% equation:
%   
resetA = ...
[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetb = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   t>=120
guard = halfspace([0,0,0,0,-1],-120);

trans{2} = transition(guard, reset, 3, 'dummy', 'names');

loc{1} = location('S1',1, inv, trans, dynamics);



%-------------------------------State P2-----------------------------------

%% equation:
%   
%               x'==vx &
%               y'==vy &
%               vx'== (n^2 + K2_11/m_c)*x + (2*n + K2_14/m_c)*vy + K2_12/m_c * y + K2_13/m_c * vx + mu/r^2 - mu/sqrt((r+x)^2-y^2)^3 * (r+x)&
%               vy'== (n^2 + K2_22/m_c)*y + (K2_23/m_c -2*n)*vx + K2_21/m_c * x + K2_24/m_c * vy - mu/sqrt((r+x)^2-y^2)^3 * y &
%               t'==1
%          
dynOpt = struct('tensorOrder',1);
dynamics = nonlinearSys(5,1,@rendeszvous_nonlinear_passive_St2_FlowEq,dynOpt); 

%% equation:
%   t<=125 & y>=-100 & x+y>=-141.1 & x>=-100 & y-x<=141.1 & y<=100 & x+y<=141.1 & x<=100 & y-x>=-141.1
invA = ...
[0,0,0,0,1;0,-1,0,0,0;-1,-1,0,0,0;-1,0,0,0,0;-1,1,0,0,0;0,1,0,0,0;1,1,0,...
0,0;1,0,0,0,0;1,-1,0,0,0];
invb = ...
[125;100;141.1;100;141.1;100;141.1;100;141.1];
invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

trans = {};
%% equation:
%   
resetA = ...
[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetb = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   t>=120
guard = halfspace([0,0,0,0,-1],-120);

trans{1} = transition(guard, reset, 3, 'dummy', 'names');

loc{2} = location('S2',2, inv, trans, dynamics);



%-----------------------------State Passive--------------------------------

%% equation:
%   
%               x'==vx & 
%               y'==vy & 
%               vx'== n^2 * x + 2*n*vy + mu/r^2 - mu/sqrt((r+x)^2-y^2)^3 * (r+x) &
%               vy'== n^2*y - 2*n*vx - mu/sqrt((r+x)^2-y^2)^3 * y &
%               t'==1
%         
dynOpt = struct('tensorOrder',1);
dynamics = nonlinearSys(5,1,@rendeszvous_nonlinear_passive_St3_FlowEq,dynOpt); 

%% equation:
%   x<=10000 & x>=-10000
invA = ...
[1,0,0,0,0;-1,0,0,0,0];
invb = ...
[10000;10000];
invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

trans = {};
loc{3} = location('S3',3, inv, trans, dynamics);



HA = hybridAutomaton(loc);


end