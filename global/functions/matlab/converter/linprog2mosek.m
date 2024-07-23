function [f,a,blc,buc,blx,bux,param,cmd] = linprog2mosek(problem)
% linprog2mosek - reformulates the problem struct from MATLAB linprog to
%    the input arguments of MOSEK msklpopt
%
% Syntax:
%    [f,a,blc,buc,blx,bux,param,cmd] = linprog2mosek(problem)
%
% Inputs:
%    problem - linear program definition, with mandatory fields
%              - problem.f (cost function min f*x)
%              - problem.Aineq (inequality constraint Aineq * x = bineq)
%              - problem.bineq (inequality constraint Aineq * x = bineq)
%              - problem.Aeq (equality constraint Aeq * x = beq)
%              - problem.beq (equality constraint Aeq * x = beq)
%              - problem.lb (lower bound for optimization variable)
%              - problem.ub (upper bound for optimization variable)
%
% Outputs:
%    f - cost function
%    a, blc, buc - two-sided constraint blc <= a*x <= buc
%    blx, bux - two-sided constraint blx <= x <= bux
%    param - empty []
%    cmd - minimization (enforced)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadprog2mosek, CORAlinprog

% Authors:       Mark Wetzlinger
% Written:       16-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only for linprog, but we don't check
%   strcmp(problem.solver,'linprog')
% for speed reasons

% note: all variables need to be here (can be empty)!
f = reshape(problem.f,[],1);
a = [problem.Aineq; problem.Aeq];
blc = [-Inf(size(problem.Aineq,1),1); problem.beq];
buc = [problem.bineq; problem.beq];
blx = problem.lb;
bux = problem.ub;

% solver parameters
param = [];
% 1. primal feasibility tolerance used by the interior-point optimizer for
%    linear problems (default: 1e-8)
% param.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-10;
% 2. relative gap termination tolerance used by the interior-point
%    optimizer for linear problems (default: 1e-8)
% param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-10;

% minimization, no display
cmd = 'minimize echo(0)';

% ------------------------------ END OF CODE ------------------------------
