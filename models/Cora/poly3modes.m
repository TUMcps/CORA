function HA = poly3modes()
% poly3modes - hybrid automaton representing the piecewise polynomial 
%   system with three discrete modes from Section 5.5 in [1] 
%
% Syntax:
%    HA = poly3modes()
%
% Inputs:
%    -
%
% Outputs:
%    HA - hybrid automaton (class hybridAutomaton)
%
% References: 
%   [1] X. Jin and et al. "Inferring Switched Nonlinear Dynamical Systems", 
%       Formal Aspects of Computing 2021
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

    % Location A ---------------------------
    
    % continuous dynamics \dot x = f(x,u)
    f = @(x,u) [-7; ...
                -x(1)];
    
    sys = nonlinearSys('locA',f);
    
    % invariant set
    inva = polytope([0 1],0);
    
    % transition 1
    guard1 = polytope([1,0],0,[0,1],0);
    reset1 = linearReset(eye(2));
    trans1 = transition(guard1,reset1,2);
    
    % transition 2
    guard2 = polytope([-1,0],0,[0,1],0);
    reset2 = linearReset(eye(2));
    trans2 = transition(guard2,reset2,3);
    
    % location object
    loc1 = location('locA',inva,[trans1;trans2],sys);
    
    
    % Location B ----------------------------
    
    % continuous dynamics \dot x = f(x,u)
    f = @(x,u) [0.5*x(1)^2 + 0.5*x(2); ...
                -9*x(1) + 3];
    
    sys = nonlinearSys('locB',f);
    
    % invariant set
    inva = polytope([0 -1; 1 0],[0;0]);
    
    % transition 1
    guard1 = polytope([1,0],0,[0,1],0);
    reset1 = linearReset(eye(2));
    trans1 = transition(guard1,reset1,1);
    
    % transition 2
    guard2 = polytope([0,-1],0,[1,0],0);
    reset2 = linearReset(eye(2));
    trans2 = transition(guard2,reset2,3);
    
    % location object
    loc2 = location('locB',inva,[trans1;trans2],sys);
    
    
    % Location C ----------------------------
    
    % continuous dynamics \dot x = f(x,u)
    f = @(x,u) [5; ...
                -0.1*x(1) - 10];
    
    sys = nonlinearSys('locC',f);
    
    % invariant set
    inva = polytope([0 -1; -1 0],[0;0]);
    
    % transition 1
    guard1 = polytope([-1,0],0,[0,1],0);
    reset1 = linearReset(eye(2));
    trans1 = transition(guard1,reset1,1);
    
    % transition 2
    guard2 = polytope([0,-1],0,[1,0],0);
    reset2 = linearReset(eye(2));
    trans2 = transition(guard2,reset2,2);
    
    % location object
    loc3 = location('locC',inva,[trans1;trans2],sys);
    
    
    % Hybrid Automaton -------------------------
    
    HA = hybridAutomaton([loc1;loc2;loc3]);

% ------------------------------ END OF CODE ------------------------------
