function HA = platoonContSX(~)


%% Generated on 20-Apr-2018

%-------------Automaton created from Component 'platoon11'-----------------

%% Interface Specification:
%   This section clarifies the meaning of state & input dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (platoon11):
%  state x := [x1; x2; x3; x4; x5; x6; x7; x8; x9]
%  input u := [u]

%--------------------------Component platoon11-----------------------------

%-------------------------------State new----------------------------------

%% equation:
%   x1' == x2 &
%   x2' == -x3 + u&
%   x3' == 1.605*x1 + 4.868*x2 -3.5754*x3 -0.8198*x4 + 0.427*x5 - 0.045*x6 -0.1942*x7 + 0.3626*x8 - 0.0946*x9 &
%   x4' == x5 &
%   x5' == x3 - x6&
%   x6' == 0.8718*x1 + 3.814*x2 -0.0754*x3 + 1.1936*x4 + 3.6258*x5  -3.2396*x6 -0.595*x7+ 0.1294*x8 -0.0796*x9 &
%   x7' == x8 &
%   x8' == x6 - x9 &
%   x9' == 0.7132*x1 + 3.573*x2 - 0.0964*x3 + 0.8472*x4 + 3.2568*x5 - 0.0876*x6 + 1.2726*x7 + 3.072*x8 - 3.1356*x9
dynA = ...
[0,1,0,0,0,0,0,0,0;0,0,-1,0,0,0,0,0,0;1.605,4.868,-3.5754,-0.8198,0.427,...
-0.045,-0.1942,0.3626,-0.0946;0,0,0,0,1,0,0,0,0;0,0,1,0,0,-1,0,0,0;...
0.8718,3.814,-0.0754,1.1936,3.6258,-3.2396,-0.595,0.1294,-0.0796;0,0,0,...
0,0,0,0,1,0;0,0,0,0,0,1,0,0,-1;0.7132,3.573,-0.0964,0.8472,3.2568,...
-0.0876,1.2726,3.072,-3.1356];
dynB = ...
[0;1;0;0;0;0;0;0;0];
dync = ...
[0;0;0;0;0;0;0;0;0];
dynamics = linearSys('linearSys', dynA, dynB, dync);

%% equation:
%   -9 <= u <= 1
invA = ...
zeros([0,9]);
invb = ...
zeros([0,1]);
invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

trans = {};
loc{1} = location('S1',1, inv, trans, dynamics);



HA = hybridAutomaton(loc);


end