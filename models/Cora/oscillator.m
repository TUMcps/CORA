function HA = oscillator()
% oscillator - hybrid automaton producing an osciallating signal taken from
%    Figure 5(a) in [1]
%
% Syntax:
%    HA = oscillator()
%
% Inputs:
%    -
%
% Outputs:
%    HA - hybrid automaton (class hybridAutomaton)
%
% References: 
%   [1] A. Gurung, M. Waga, and K. Suenaga. "Learning nonlinear hybrid 
%       automata from input–output time-series data", ATVA 2023
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton

% Authors:       Niklas Kochdumper
% Written:       04-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % Location 1 -------------------------------

    % system dynamics \dot x = A*x + B*u + c
    A = [-2 0;0 -1];
    c = [1.4;-0.7];
    
    linSys = linearSys(A,[],c);
    
    % invariant set
    inva = polytope([-1 0; -0.714286 -1],[0;0]);
    
    % transition
    guard = polytope([],[],[-0.714286 -1],0);
    
    reset = linearReset(eye(2));
    
    trans = transition(guard,reset,2);
    
    % location object
    loc1 = location(inva,trans,linSys);
    
    
    % Location 2 ---------------------------------
    
    % system dynamics \dot x = A*x + B*u + c
    A = [-2 0;0 -1];
    c = [-1.4;0.7];
    
    linSys = linearSys(A,[],c);
    
    % invariant set
    inva = polytope([-1 0; 0.714286 1],[0;0]);
    
    % transition
    guard = polytope([],[],[1 0],0);
    
    reset = linearReset(eye(2));
    
    trans = transition(guard,reset,3);
    
    % location object
    loc2 = location(inva,trans,linSys);
    
    
    % Location 3 ----------------------------------
    
    % system dynamics \dot x = A*x + B*u + c
    A = [-2 0;0 -1];
    c = [-1.4;0.7];
    
    linSys = linearSys(A,[],c);
    
    % invariant set
    inva = polytope([1 0; 0.714286 1],[0;0]);
    
    % transition
    guard = polytope([],[],[0.714286 1],0);
    
    reset = linearReset(eye(2));
    
    trans = transition(guard,reset,4);
    
    % location object
    loc3 = location(inva,trans,linSys);
    
    
    % Location 4 ----------------------------------
    
    % system dynamics \dot x = A*x + B*u + c
    A = [-2 0;0 -1];
    c = [1.4;-0.7];
    
    linSys = linearSys(A,[],c);
    
    % invariant set
    inva = polytope([1 0; -0.714286 -1],[0;0]);
    
    % transition
    guard = polytope([],[],[1 0],0);
    
    reset = linearReset(eye(2));
    
    trans = transition(guard,reset,1);
    
    % location object
    loc4 = location(inva,trans,linSys);
    
    
    % Hybrid Automaton ---------------------------
    
    HA = hybridAutomaton([loc1;loc2;loc3;loc4]);

% ------------------------------ END OF CODE ------------------------------
