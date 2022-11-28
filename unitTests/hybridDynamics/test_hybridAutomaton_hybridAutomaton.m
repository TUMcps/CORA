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

% Author:       Mark Wetzlinger
% Written:      26-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% invariant
polyOpt = struct('A',[-1,0],'b',0);
inv_2D = mptPolytope(polyOpt);
polyOpt = struct('A',[-1,0,0],'b',0);
inv_3D = mptPolytope(polyOpt);

% transition
c = [-1;0]; d = 0; C = [0,1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
trans_2D{1} = transition(guard,reset,2);
c = [-1;0;0]; d = 0; C = [0,0,1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',eye(3),'c',zeros(3,1));
trans_3D{1} = transition(guard,reset,1);

% flow equation
dynamics_2D = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics_3D = linearSys([0,1,0;0,0,1;0,0,1],1);

% define location
loc_1 = location('S1',inv_2D,trans_2D,dynamics_2D);
loc_2 = location('S2',inv_3D,trans_3D,dynamics_3D);


try
    % dimensions of reset functions do not match
    loc = {loc_1,loc_2};
    HA = hybridAutomaton(loc);
    res = false;
end
try
    % too many input arguments
    loc{1} = loc_1;
    HA = hybridAutomaton(loc,loc);
    res = false;
end
try
    % not a cell array
    HA = hybridAutomaton(1);
    res = false;
end
try
    % not a cell array
    HA = hybridAutomaton(loc_1);
    res = false;
end
try
    % not a cell array of locations
    HA = hybridAutomaton({'wrong',dynamics_2D});
    res = false;
end
try
    % not a cell array of all locations
    HA = hybridAutomaton([{loc_1}; {1}]);
    res = false;
end

%------------- END OF CODE --------------
