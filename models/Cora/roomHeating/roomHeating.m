function HA = roomHeating()
% roomHeating - room heating benchmark described in Sec. 2.3 in [1] with 
%               two rooms 
%
% Syntax:  
%    PHA = roomHeating()
%
% Inputs:
%    ---
%
% Outputs:
%    HA - hybridAutomaton object
%
% References:
%   [1] A. Fehnker and F. Ivancic. "Benchmarks for Hybrid Systems 
%       Verification", HSCC 2004

% Author:       Niklas Kochdumper
% Written:      26-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Test parallelHybridAutomaton

% parameter room 1
a1 = 0.5; 
b1 = 0.4; 
c1 = 6;

% parameter room 2
a2 = 0.5;
b2 = 0.3;
c2 = 7;

% temperatures for switching heaters on and off
T_off = 21;
T_on = 20;


% Location 1 - room 1 on, room 2 on ---------------------------------------

% system dynamics
A = [-(a1 + b1), a1; a2 -(a2 + b2)];
B = [b1;b2];
c = [c1;c2];

linSys = linearSys(A,B,c);

% invariant set
inv = polytope([1 0; 0 1],[T_off;T_off]);

% transition 1: room 1 on -> off
guard = conHyperplane([1 0],T_off);

reset.A = eye(2);
reset.c = zeros(2,1);

trans = transition(guard,reset,2);

% transition 2: room 2 on -> off
guard = conHyperplane([0 1],T_off);

reset.A = eye(2);
reset.c = zeros(2,1);

trans(2) = transition(guard,reset,3);

% location object
loc = location('on',inv,trans,linSys);


% Location 2 - room 1 off, room 2 on --------------------------------------

% system dynamics
A = [-(a1 + b1), a1; a2 -(a2 + b2)];
B = [b1;b2];
c = [0;c2];

linSys = linearSys(A,B,c);

% invariant set
inv = polytope([-1 0; 0 1],[-T_on;T_off]);

% transition 1: room 1 off -> on
guard = conHyperplane([1 0],T_on);

reset.A = eye(2);
reset.c = zeros(2,1);

trans(1) = transition(guard,reset,1);

% transition 2: room 2 on -> off
guard = conHyperplane([0 1],T_off);

reset.A = eye(2);
reset.c = zeros(2,1);

trans(2) = transition(guard,reset,4);

% location object
loc(2) = location('on',inv,trans,linSys);


% Location 3 - room 1 on, room 2 off --------------------------------------

% system dynamics
A = [-(a1 + b1), a1; a2 -(a2 + b2)];
B = [b1;b2];
c = [c1;0];

linSys = linearSys(A,B,c);

% invariant set
inv = polytope([1 0; 0 -1],[T_off;-T_on]);

% transition 1: room 1 on -> off
guard = conHyperplane([1 0],T_off);

reset.A = eye(2);
reset.c = zeros(2,1);

trans(1) = transition(guard,reset,4);

% transition 2: room 2 off -> on
guard = conHyperplane([0 1],T_on);

reset.A = eye(2);
reset.c = zeros(2,1);

trans(2) = transition(guard,reset,1);

% location object
loc(3) = location('on',inv,trans,linSys);


% Location 4 - room 1 off, room 2 off -------------------------------------

% system dynamics
A = [-(a1 + b1), a1; a2 -(a2 + b2)];
B = [b1;b2];

linSys = linearSys(A,B);

% invariant set
inv = polytope([-1 0; 0 -1],[-T_on;-T_on]);

% transition 1: room 1 on -> off
guard = conHyperplane([1 0],T_on);

reset.A = eye(2);
reset.c = zeros(2,1);

trans(1) = transition(guard,reset,2);

% transition 2: room 2 on -> off
guard = conHyperplane([0 1],T_on);

reset.A = eye(2);
reset.c = zeros(2,1);

trans(2) = transition(guard,reset,3);

% location object
loc(4) = location('on',inv,trans,linSys);


% Hybrid Automaton --------------------------------------------------------

HA = hybridAutomaton(loc);

end

%------------- END OF CODE --------------