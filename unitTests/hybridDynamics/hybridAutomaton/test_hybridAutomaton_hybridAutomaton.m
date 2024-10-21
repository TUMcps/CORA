function res = test_hybridAutomaton_hybridAutomaton
% test_hybridAutomaton_hybridAutomaton - unit test function for constructor
%
% Syntax:
%    res = test_hybridAutomaton_hybridAutomaton
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Mark Wetzlinger
% Written:       26-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty object
hybridAutomaton()

% invariant
inv_2D = polytope([-1,0],0);
inv_3D = polytope([-1,0,0],0);

% transition
guard = polytope([0 1],0,[-1 0],0);
reset = linearReset([1,0;0,-0.75]);
trans_2D = transition(guard,reset,2);
guard = polytope([0,0,1],0,[-1 0 0],0);
reset = linearReset.eye(3);
trans_3D = transition(guard,reset,1);

% flow equation
dynamics_2D = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics_3D = linearSys([0,1,0;0,0,1;0,0,1],1);

% define location
loc_1 = location('S1',inv_2D,trans_2D,dynamics_2D);
loc_2 = location('S2',inv_3D,trans_3D,dynamics_3D);

% init hybrid automaton
HA = hybridAutomaton([loc_1;loc_1]);
% init with name
HA = hybridAutomaton('HA',[loc_1;loc_1]);


% wrong initializations

% wrong name
assertThrowsAs(@hybridAutomaton,'CORA:wrongValue',2,[loc_1;loc_1]);

% dimensions of reset functions do not match
assertThrowsAs(@hybridAutomaton,'CORA:wrongInputInConstructor',[loc_1;loc_2]);

% too many input arguments
assertThrowsAs(@hybridAutomaton,'CORA:numInputArgsConstructor','HA',loc_1,loc_1);

% not a location array
assertThrowsAs(@hybridAutomaton,'CORA:wrongInputInConstructor',{loc_1,loc_2});


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
