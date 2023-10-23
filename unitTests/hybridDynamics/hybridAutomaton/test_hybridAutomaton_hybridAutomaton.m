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

% assume true
res = true;

% empty object
hybridAutomaton()

% invariant
inv_2D = polytope([-1,0],0);
inv_3D = polytope([-1,0,0],0);

% transition
c = [-1;0]; d = 0; C = [0,1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
trans_2D = transition(guard,reset,2);
c = [-1;0;0]; d = 0; C = [0,0,1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',eye(3),'c',zeros(3,1));
trans_3D = transition(guard,reset,1);

% flow equation
dynamics_2D = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics_3D = linearSys([0,1,0;0,0,1;0,0,1],1);

% define location
loc_1 = location('S1',inv_2D,trans_2D,dynamics_2D);
loc_2 = location('S2',inv_3D,trans_3D,dynamics_3D);

% init hybrid automaton
HA = hybridAutomaton([loc_1;loc_1]);


% wrong initializations
if CHECKS_ENABLED

try
    % dimensions of reset functions do not match
    loc = [loc_1;loc_2];
    HA = hybridAutomaton(loc);
    res = false;
end
try
    % too many input arguments
    HA = hybridAutomaton(loc_1,loc_1);
    res = false;
end
try
    % not a location array
    HA = hybridAutomaton({loc_1,loc_2});
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
