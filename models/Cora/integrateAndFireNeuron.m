function HA = integrateAndFireNeuron()
% integrateAndFireNeuron - integrate and fire neuron modelled as a hybrid
%                          automaton with the three locations "charging",
%                          "firing", and "reset"
%
% Syntax:  
%    HA = integrateAndFireNeuron()
%
% Inputs:
%    no
%
% Outputs:
%    HA - object of class hybridAutomaton

% Author:       Niklas Kochdumper
% Written:      15-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


    % Parameter -----------------------------------------------------------
    
    % the system state is [u;t], where u is the voltage and t is the time
    
    Cm = 1e-6;                  % membrane capacity
    R = 1e4;                    % membrane resistance
    u_rest = 0;                 % equilibrium voltage
    u_fire = 40e-3;             % firing voltage threshold
    u_out = 80e-3;              % output voltage
    
    
    % Location 1: Charging ------------------------------------------------
    
    % system dynamics
    A = [-1/(R*Cm) 0; 0 0];
    B = [1/Cm 1/Cm; 0 0];
    c = [u_rest/(R*Cm); 1];
    C = [1/R 0];
    
    sys = linearSys(A,B,c,C);
    
    % invariant set
    inv = polytope([1 0],u_fire);
    
    % transition
    guard = conHyperplane([1 0],u_fire);
    reset.A = eye(2); reset.c = [0;0];
    
    tran = transition(guard,reset,2);
    
    % location object
    loc(1) = location(inv,tran,sys);
    
    
    % Location 2: Firing --------------------------------------------------

    % system dynamics
    A = [-1/(R*Cm) 0; 0 0];
    B = [1/Cm 1/Cm; 0 0];
    c = [u_fire/(R*Cm) + 1e-4/Cm; 1];
    C = [1/R 0];
    
    sys = linearSys(A,B,c,C);
    
    % invariant set
    inv = polytope([1 0],u_out);
    
    % transition
    guard = conHyperplane([1 0],u_out);
    reset.A = eye(2); reset.c = [0;0];
    
    tran = transition(guard,reset,3);
    
    % location object
    loc(2) = location(inv,tran,sys);
    
    
    % Location 3: Reset ---------------------------------------------------

    % system dynamics
    A = [-1/(R*Cm) 0; 0 0];
    B = [1/Cm 1/Cm; 0 0];
    c = [u_rest/(R*Cm) - 1e-4/Cm; 1];
    C = [1/R 0];
    
    sys = linearSys(A,B,c,C);
    
    % invariant set
    inv = polytope([-1 0],-u_rest);
    
    % transition
    guard = conHyperplane([1 0],u_rest);
    reset.A = eye(2); reset.c = [0;0];
    
    tran = transition(guard,reset,1);
    
    % location object
    loc(3) = location(inv,tran,sys);

    
    % Hybrid Automaton ----------------------------------------------------
    
    HA = hybridAutomaton(loc);

end

%------------- END OF CODE --------------