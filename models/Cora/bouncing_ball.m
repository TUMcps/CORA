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
%    HA - hybrid automaton

% Author:       ---
% Written:      ---
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% continuous dynamics
A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];
linSys = linearSys('linearSys',A,B,c);

% invariant set
inv = polytope([-1,0],0);

% guard sets
guard = conHyperplane([1,0],0,[0,1],0);

% reset function
reset.A = [0, 0; 0, alpha];
reset.c = zeros(2,1);

% transitions
trans = transition(guard,reset,1);

% location object
loc = location('loc1',inv,trans,linSys);

% hybrid automata
HA = hybridAutomaton(loc);

end

%------------- END OF CODE --------------
