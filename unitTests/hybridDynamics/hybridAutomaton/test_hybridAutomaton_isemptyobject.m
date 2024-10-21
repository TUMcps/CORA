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
assert(isemptyobject(hybridAutomaton()));

% invariant
inv_2D_poly = polytope([-1,0],0);
inv_3D_int = interval([-1;0;0],[1;2;3]);

% transitions
guard = polytope([0 1],0,[-1 0],0);
reset = linearReset([1,0;0,-0.75;0,-1]);
trans_2D_3D = transition(guard,reset,2);
guard = polytope([0 0 1],0,[-1 0 0],0);
reset = linearReset(zeros(2,3));
trans_3D_2D = transition(guard,reset,1);

% flow equation
dynamics_2D_lin = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics_3D_lin = linearSys([0,1,0;0,0,1;0,0,1],1);

% define locations
loc1 = location(inv_2D_poly,trans_2D_3D,dynamics_2D_lin);
loc2 = location(inv_3D_int,trans_3D_2D,dynamics_3D_lin);

% non-empty hybrid automaton
HA = hybridAutomaton([loc1;loc2]);
assert(~isemptyobject(HA));

% array of hybrid automata
assert(all(isemptyobject([hybridAutomaton(),HA]) == [true false]));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
