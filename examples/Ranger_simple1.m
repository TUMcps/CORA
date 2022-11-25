function HA = Ranger_simple1(~)


%% Generated on 24-May-2022

%--------------Automaton created from Component 'Movement'-----------------

%% Interface Specification:
%   This section clarifies the meaning of state, input & output dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (Movement):
%  state x := [pos; vel; acc; t]
%  input u := [const_acc; bound]

%--------------------------Component Movement------------------------------

%----------------------------State Movement--------------------------------

%% equation:
%   pos' == vel &
%   vel' == acc &
%   acc' == 0 &
%   t' == 1
dynA = ...
[0,1,0,0;0,0,1,0;0,0,0,0;0,0,0,0];
dynB = ...
[0,0;0,0;0,0;0,0];
dync = ...
[0;0;0;1];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   pos <= bound &
%   pos >= -bound
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   acc := -const_acc &
%   pos := bound &
%   vel := 0
resetA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,1];
resetB = ...
[0,1;0,0;-1,0;0,0];
resetc = ...
[0;0;0;0];
resetInputDim = 2;
reset = struct('A', resetA, 'B', resetB,'c', resetc,'hasInput',1,'inputDim',resetInputDim);

%% equation:
%   pos >= bound
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{1} = transition(guard, reset, 1,[], 4);

%% equation:
%   acc := const_acc &
%   pos := -bound &
%   vel := 0
resetA = ...
[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,1];
resetB = ...
[0,-1;0,0;1,0;0,0];
resetc = ...
[0;0;0;0];
resetInputDim = 2;
reset = struct('A', resetA, 'B', resetB,'c', resetc,'hasInput',1,'inputDim',resetInputDim);

%% equation:
%   pos <= -bound
A = ...
zeros([0,0]);
b = ...
zeros([0,1]);
polyOpt = struct('A', A, 'b', b);
guard = mptPolytope(polyOpt);

trans{2} = transition(guard, reset, 1,[], 4);

loc{1} = location('S1', inv, trans, dynamics);



HA = hybridAutomaton(loc);


end