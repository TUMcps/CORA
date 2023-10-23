function res = test_hybridAutomaton_isemptyobject
% test_hybridAutomaton_isemptyobject - test function for emptiness check
%
% Syntax:
%    res = test_hybridAutomaton_isemptyobject
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty automaton
res = isemptyobject(hybridAutomaton());

% invariant
inv_2D_poly = polytope([-1,0],0);
inv_3D_int = interval([-1;0;0],[1;2;3]);

% transitions
c = [-1;0]; d = 0; C = [0,1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',[1,0;0,-0.75;0,-1],'c',[0;0;0]);
trans_2D_3D = transition(guard,reset,2);
c = [-1;0;0]; d = 0; C = [0,0,1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',zeros(2,3),'c',zeros(2,1));
trans_3D_2D = transition(guard,reset,1);

% flow equation
dynamics_2D_lin = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics_3D_lin = linearSys([0,1,0;0,0,1;0,0,1],1);

% define locations
loc1 = location(inv_2D_poly,trans_2D_3D,dynamics_2D_lin);
loc2 = location(inv_3D_int,trans_3D_2D,dynamics_3D_lin);

% non-empty hybrid automaton
HA = hybridAutomaton([loc1;loc2]);
res(end+1,1) = ~isemptyobject(HA);

% array of hybrid automata
res(end+1,1) = all(isemptyobject([hybridAutomaton(),HA]) == [true false]);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
