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
%    [1] L. Geretti, “ARCH-COMP20 Category Report: Continuous and Hybrid 
%        Systems with Nonlinear Dynamics", 2020

% Author:        Niklas Kochdumper
% Written:       19-June-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % Location Outside ----------------------------------------------------

    % dynamics
    sys = nonlinearSys(@lotkaVolterraOutside);

    % invariant set
    syms x y c t
    vars = [x;y;c;t];
    eq = -(x-1)^2 - (y-1)^2 + 0.15^2;

    inv = levelSet(eq,vars,'<=');

    % transition
    guard = levelSet(eq,vars,'==');
    reset.A = eye(4); reset.b = zeros(4,1);

    trans{1} = transition(guard, reset, 2);

    % location
    loc{1} = location('outside', inv, trans, sys);


    % Location Inside -----------------------------------------------------

    % dynamics
    sys = nonlinearSys(@lotkaVolterraInside);

    % invariant set
    syms x y c t
    vars = [x;y;c;t];
    eq = (x-1)^2 + (y-1)^2 - 0.15^2;

    inv = levelSet(eq,vars,'<=');

    % location
    loc{2} = location('inside', inv, [], sys);


    % Hybrid Automaton ----------------------------------------------------

    HA = hybridAutomaton(loc);
    
end

% Auxiliary Functions -----------------------------------------------------

function f = lotkaVolterraInside(x,u)

    f(1,1) = 3*x(1) - 3*x(1)*x(2);
    f(2,1) = x(1)*x(2) - x(2);
    f(3,1) = 1;
    f(4,1) = 1;
end

function f = lotkaVolterraOutside(x,u)

    f(1,1) = 3*x(1) - 3*x(1)*x(2);
    f(2,1) = x(1)*x(2) - x(2);
    f(3,1) = 0;
    f(4,1) = 1;
end

%------------- END OF CODE --------------