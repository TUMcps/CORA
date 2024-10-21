function example_manual_bouncing_ball_construction()
% example_manual_bouncing_ball_construction - example from the manual 
%    demonstrating the construction of the bouncing ball example
%
% Syntax:
%    example_manual_bouncing_ball_construction()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% guard set
guard = polytope([0 1],0,[1 0],0);

% reset function
reset = linearReset([1 0; 0 -0.75],[0;0],[0;0]);

% transtition object
trans = transition(guard,reset,1);

% differential equation
sys = linearSys([0 1;0 0],[0;0],[0;-9.81]);

% invariant set
inv = polytope([-1 0],0);

% location object
loc = location(inv,trans,sys);

% list of locations
locs(1) = loc;

% hybrid automaton object
HA = hybridAutomaton('bouncingBall',locs);

% ------------------------------ END OF CODE ------------------------------
