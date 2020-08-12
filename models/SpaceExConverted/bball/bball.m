function HA = bball(~)


%% Generated on 15-Jul-2020

%---------------Automaton created from Component 'system'------------------

%% Interface Specification:
%   This section clarifies the meaning of state, input & output dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (system.ball):
%  state x := [x; v]
%  input u := [uDummy]

%-------------------------Component system.ball----------------------------

%-----------------------------State always---------------------------------

%% equation:
%   x' == v & v' == -g
dynA = ...
[0,1;0,0];
dynB = ...
[0;0];
dync = ...
[0;-9.81];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   x >= 0
A = ...
[-1,0];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   v' := -c*v
resetA = ...
[1,0;0,-0.75];
resetb = ...
[0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   x <= eps & v < 0
c = [-1;0];
d = 0;C = ...
[0,1];
D = [0];

guard = conHyperplane(c,d,C,D);

trans{1} = transition(guard, reset, 1);

loc{1} = location('S1', inv, trans, dynamics);



HA = hybridAutomaton(loc);


end