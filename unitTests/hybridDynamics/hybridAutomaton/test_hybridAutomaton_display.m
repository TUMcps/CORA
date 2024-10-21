function res = test_hybridAutomaton_display
% test_hybridAutomaton_display - test function for display
%
% Syntax:
%    res = test_hybridAutomaton_display
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
% Written:       08-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty object
hybridAutomaton()

% name
name = 'HA';

% single location
inv = polytope([-1,0],0);
guard = polytope([0 1],0,[-1 0],0);
reset = linearReset([1,0;0,-0.75]);
trans(1) = transition(guard,reset,1);
dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
loc(1) = location('S1',inv,trans,dynamics);
HA = hybridAutomaton(name,loc)


% multiple locations
inv = polytope([-1,0,0],0);
dynamics = linearSys([1 -2 1; 0 0 2; -1 0 1]);
guard = polytope([0,0,1],0,[-1 0 0],0);
reset = linearReset([1 0 0; 0 1 0]);
trans(1) = transition(guard,reset,2);
guard = polytope([1 0 0],0,[-1 0 0],0);
reset = linearReset([-1 0 0; 0 1 0]);
trans(2) = transition(guard,reset,2);
loc(1) = location('S1',inv,trans,dynamics);

inv = polytope([-1,0],0);
dynamics = linearSys([-2 1; 0 -1]);
guard = polytope([0 1],0,[-1 0],0);
reset = linearReset([1 0;0 1;0 2],[],[-1;0;0]);
trans_(1) = transition(guard,reset,1);
loc(2) = location('S1',inv,trans_,dynamics);

HA = hybridAutomaton(name,loc)


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
