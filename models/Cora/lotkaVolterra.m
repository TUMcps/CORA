function HA = lotkaVolterra()
% lotkaVolterra - hybrid Lotka-Volterra model with tangential guard crossing
%               (see Sec. 3.5 in [1])
%
% Syntax:  
%    HA = lotkaVolterra()
%
% Inputs:
%    ---
%
% Outputs:
%    HA - hybridAutomaton object
% 
% References:
%    [1] L. Geretti et al., "ARCH-COMP20 Category Report: Continuous and
%        Hybrid Systems with Nonlinear Dynamics", 2020

% Author:        Niklas Kochdumper
% Written:       19-June-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % Location Outside ----------------------------------------------------
    
    % dynamics
    sys = nonlinearSys(@lotkaVolterraDyn);

    % invariant set
    syms x y t
    vars = [x;y;t];
    eq = -(x-1)^2 - (y-1)^2 + 0.161^2;

    inv = levelSet(eq,vars,'<=');

    % transition
    guard = levelSet(eq,vars,'==');
    reset.A = eye(3); reset.c = zeros(3,1);

    trans = transition(guard, reset, 2);

    % location
    loc(1) = location('outside', inv, trans, sys);


    % Location Inside -----------------------------------------------------

    % dynamics
    sys = nonlinearSys(@lotkaVolterraDyn);

    % invariant set
    syms x y t
    vars = [x;y;t];
    eq = (x-1)^2 + (y-1)^2 - 0.161^2;

    inv = levelSet(eq,vars,'<=');

    % location
    loc(2) = location('inside', inv, transition(), sys);
    
    % cannot model transition because of infinite switching
    % transition
%     guard = levelSet(eq,vars,'==');
%     reset.A = eye(3); reset.c = zeros(3,1);
% 
%     trans = transition(guard, reset, 1);
% 
%     % location
%     loc(2) = location('inside', inv, trans, sys);


    % Hybrid Automaton ----------------------------------------------------

    HA = hybridAutomaton(loc);
    
end

% Auxiliary Functions -----------------------------------------------------

function f = lotkaVolterraDyn(x,u)

    f = [3*x(1) - 3*x(1)*x(2);
         x(1)*x(2) - x(2);
         1];
    
end

%------------- END OF CODE --------------