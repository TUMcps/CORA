function HA = lotkaVolterra()
% lotkaVolterra - hybrid Lotka-Volterra model with tangential guard
%    crossing (see Sec. 3.5 in [1])
%
% Syntax:
%    HA = lotkaVolterra()
%
% Inputs:
%    -
%
% Outputs:
%    HA - hybridAutomaton object
% 
% References:
%    [1] L. Geretti et al., "ARCH-COMP20 Category Report: Continuous and
%        Hybrid Systems with Nonlinear Dynamics", 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       19-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Location Outside --------------------------------------------------------

% dynamics
sys = nonlinearSys(@aux_lotkaVolterraDyn);

% invariant set
syms x y t
vars = [x;y;t];
eq = -(x-1)^2 - (y-1)^2 + 0.161^2;

inv = levelSet(eq,vars,'<=');

% transition
guard = levelSet(eq,vars,'==');
reset = linearReset.eye(3);

trans = transition(guard, reset, 2);

% location
loc(1) = location('outside', inv, trans, sys);


% Location Inside ---------------------------------------------------------

% dynamics
sys = nonlinearSys(@aux_lotkaVolterraDyn);

% invariant set
syms x y t
vars = [x;y;t];
eq = (x-1)^2 + (y-1)^2 - 0.161^2;

inv = levelSet(eq,vars,'<=');

% location
loc(2) = location('inside', inv, transition(), sys);

% cannot model transition because of infinite switching transition
% guard = levelSet(eq,vars,'==');
% reset = linearReset(eye(3),zeros(3,1),zeros(3,1));
% 
% trans = transition(guard, reset, 1);
% 
% % location
% loc(2) = location('inside', inv, trans, sys);


% compose hybrid automaton
HA = hybridAutomaton('lotkaVolterra',loc);
    
end


% Auxiliary functions -----------------------------------------------------

function f = aux_lotkaVolterraDyn(x,u)

    f = [3*x(1) - 3*x(1)*x(2);
         x(1)*x(2) - x(2);
         1];
    
end

% ------------------------------ END OF CODE ------------------------------
