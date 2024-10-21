function HA = bouncing_ball(alpha)
% bouncing_ball - bouncing ball system from the manual
%
% Syntax:
%    HA = bouncing_ball(alpha)
%
% Inputs:
%    alpha - damping factor
%
% Outputs:
%    HA - hybridAutomaton object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ---
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% continuous dynamics
A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];
linsys = linearSys('linearSys',A,B,c);

% invariant set
inv = polytope([-1,0],0);

% guard sets
guard = polytope([0,1],0,[1,0],0);

% reset function
reset = linearReset([0, 0; 0, alpha],zeros(2,1),zeros(2,1));

% transitions
trans = transition(guard,reset,1);

% location object
loc = location('loc1',inv,trans,linsys);

% hybrid automata
HA = hybridAutomaton('BouncingBall',loc);

% ------------------------------ END OF CODE ------------------------------
