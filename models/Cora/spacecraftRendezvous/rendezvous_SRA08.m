function HA = rendezvous_SRA08(~)
% rendezvous_SRA08 - linear spacecraft-rendezvous benchmark 
%                    (see Sec. 3.2 in [1])
%
% Syntax:  
%    HA = rendezvous_SRA08()
%
% Inputs:
%    ---
%
% Outputs:
%    HA - hybridAutomaton object
% 
% References:
%    [1] M. Althoff, “ARCH-COMP19 Category Report: Continuous and Hybrid 
%        Systems with Linear Continuous Dynamics", 2019

% Author:        Niklas Kochdumper
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

%% Generated on 20-Apr-2018

%----------Automaton created from Component 'ChaserSpacecraft'-------------

%% Interface Specification:
%   This section clarifies the meaning of state & input dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (ChaserSpacecraft):
%  state x := [x; y; vx; vy; t]
%  input u := [uDummy]

%----------------------Component ChaserSpacecraft--------------------------

% times for transition to aborting mode
t_l = 0;
t_u = 240;
dist = 1e5;

%-------------------------------State P2-----------------------------------

%% equation:
%   t'==1 & x'==vx & y'==vy & vx'==-0.057599765881773*x+0.000200959896519766*y-2.89995083970656*vx+0.00877200894463775*vy & vy'==-0.000174031357370456*x-0.0665123984901026*y-0.00875351105536225*vx-2.90300269286856*vy
dynA = ...
[0,0,1,0,0;0,0,0,1,0;-0.05759976588,0.0002009598965,-2.89995084,...
0.008772008945,0;-0.0001740313574,-0.06651239849,-0.008753511055,...
-2.903002693,0;0,0,0,0,0];
dynB = ...
[0;0;0;0;0];
dync = ...
[0;0;0;0;1];
dynamics = linearSys('linearSys', dynA, dynB, dync);

%% equation:
%   t<=t_u & x<=-100
invA = ...
[0,0,0,0,1;1,0,0,0,0];
invb = ...
[t_u;-100];
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
%   y>=-100 & x+y >=-141.1 & x>=-100 & y-x<=141.1 & y<=100 & x+y<=141.1 & x<=100 & y-x>=-141.1
guardA = [-1;0;0;0;0];
guardb = 100;
guard = conHyperplane(guardA,guardb);

trans(1) = transition(guard, reset, 2);

%% equation:
%   
resetA = ...
[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetc = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'c', resetc);

%% equation:
%   t>=t_l & t<=0.1
guard = interval([-dist;-dist;-dist;-dist;t_l],[dist;dist;dist;dist;0.1]);
trans(2) = transition(guard, reset, 3);

%   t>=0.1 & t<=0.15
guard = interval([-dist;-dist;-dist;-dist;0.1],[dist;dist;dist;dist;0.15]);
trans(3) = transition(guard, reset, 3);

%   t>=0.15 & t<=0.2
guard = interval([-dist;-dist;-dist;-dist;0.15],[dist;dist;dist;dist;0.2]);
trans(4) = transition(guard, reset, 3);

%   t>=0.2 & t<=3
guard = interval([-dist;-dist;-dist;-dist;0.2],[dist;dist;dist;dist;3]);
trans(5) = transition(guard, reset, 3);

%   t>=3 & t<=t_u
guard = interval([-dist;-dist;-dist;-dist;3],[dist;dist;dist;dist;t_u]);
trans(6) = transition(guard, reset, 3);

loc(1) = location('S1', inv, trans, dynamics);



%-------------------------------State P3-----------------------------------

%% equation:
%   t'==1 & x'==vx & y'==vy & vx'==-0.575999943070835*x+0.000262486079431672*y-19.2299795908647*vx+0.00876275931760007*vy & vy'==-0.000262486080737868*x-0.575999940191886*y-0.00876276068239993*vx-19.2299765959399*vy
dynA = ...
[0,0,1,0,0;0,0,0,1,0;-0.5759999431,0.0002624860794,-19.22997959,...
0.008762759318,0;-0.0002624860807,-0.5759999402,-0.008762760682,...
-19.2299766,0;0,0,0,0,0];
dynB = ...
[0;0;0;0;0];
dync = ...
[0;0;0;0;1];
dynamics = linearSys('linearSys', dynA, dynB, dync);

%% equation:
%   t<=t_u & y>=-100 & x+y>=-141.1 & x>=-100 & y-x<=141.1 & y<=100 & x+y<=141.1 & x<=100 & y-x>=-141.1
invA = ...
[0,0,0,0,1;0,-1,0,0,0;-1,-1,0,0,0;-1,0,0,0,0;-1,1,0,0,0;0,1,0,0,0;1,1,0,...
0,0;1,0,0,0,0;1,-1,0,0,0];
invb = ...
[t_u;100;141.1;100;141.1;100;141.1;100;141.1];
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
%  t>=t_l & t<=110
guard = interval([-dist;-dist;-dist;-dist;t_l],[dist;dist;dist;dist;110]);
trans(1) = transition(guard, reset, 3);

%  t>=110 & t<=112
guard = interval([-dist;-dist;-dist;-dist;110],[dist;dist;dist;dist;112]);
trans(2) = transition(guard, reset, 3);

%  t>=112 & t<=115
guard = interval([-dist;-dist;-dist;-dist;112],[dist;dist;dist;dist;115]);
trans(3) = transition(guard, reset, 3);

%  t>=115 & t<=120
guard = interval([-dist;-dist;-dist;-dist;115],[dist;dist;dist;dist;120]);
trans(4) = transition(guard, reset, 3);

%  t>=120 & t<=135
guard = interval([-dist;-dist;-dist;-dist;120],[dist;dist;dist;dist;135]);
trans(5) = transition(guard, reset, 3);

%  t>=135 & t<=t_u
guard = interval([-dist;-dist;-dist;-dist;135],[dist;dist;dist;dist;t_u]);
trans(6) = transition(guard, reset, 3);

loc(2) = location('S2',inv, trans, dynamics);



%-----------------------------State Passive--------------------------------

%% equation:
%   t'==1 & x'==vx & y'==vy & vx'==0.0000575894721132000*x+0.00876276*vy & vy'==-0.00876276*vx 
dynA = ...
[0,0,1,0,0;0,0,0,1,0;5.758947211E-05,0,0,0.00876276,0;0,0,-0.00876276,0,...
0;0,0,0,0,0];
dynB = ...
[0;0;0;0;0];
dync = ...
[0;0;0;0;1];
dynamics = linearSys('linearSys', dynA, dynB, dync);

%% equation:
%   t>=40
invA = ...
[0,0,0,0,-1];
invb = ...
40;
inv = polytope(invA, invb);

trans = transition();
loc(3) = location('S3',inv, trans, dynamics);



HA = hybridAutomaton(loc);


end

%------------- END OF CODE --------------