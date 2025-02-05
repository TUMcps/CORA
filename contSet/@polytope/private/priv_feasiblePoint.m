function x = priv_feasiblePoint(P)
% priv_feasiblePoint - computes a feasible point of the polytope P, if it
%   exists
%
% Syntax:
%    x = priv_feasiblePoint(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    x - feasible point x\in P. If no such point exists, return x = []
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       17-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init linprog struct
problem.f = zeros([dim(P) 1]); % Dummy function
problem.Aineq = P.A_.val;
problem.bineq = P.b_.val;
problem.Aeq = P.Ae_.val;
problem.beq = P.be_.val;
problem.lb = [];
problem.ub = [];

x = CORAlinprog(problem);

% ------------------------------ END OF CODE ------------------------------
