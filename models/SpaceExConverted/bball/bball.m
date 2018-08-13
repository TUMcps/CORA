function HA = bball(~)


%% Generated on 09-Aug-2018

%---------------Automaton created from Component 'system'------------------

%% Interface Specification:
%   This section clarifies the meaning of state & input dimensions
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
[0;-9.8100000000000004973799150320701];
dynamics = linearSys('linearSys', dynA, dynB, dync);

%% equation:
%   x >= 0
invA = ...
[-1,0];
invb = ...
[-0];
invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

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
guardA = ...
[1,0;0,1];
guardb = ...
[-0;-0];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{1} = transition(guard, reset, 1, 'dummy', 'names');

loc{1} = location('S1',1, inv, trans, dynamics);



HA = hybridAutomaton(loc);


end