function [q,c,a,blc,buc,blx,bux,param,cmd] = quadprog2mosek(problem)
% quadprog2mosek - reformulates the problem struct from MATLAB quadprog to
%    the input arguments of MOSEK msklpopt
%
% Syntax:
%    [q,c,a,blc,buc,blx,bux,param,cmd] = quadprog2mosek(problem)
%
% Inputs:
%    problem - quadratic program definition, with mandatory fields
%              - problem.f (cost function min f*x)
%              - problem.Aineq (inequality constraint Aineq * x = bineq)
%              - problem.bineq (inequality constraint Aineq * x = bineq)
%              - problem.Aeq (equality constraint Aeq * x = beq)
%              - problem.beq (equality constraint Aeq * x = beq)
%              - problem.lb (lower bound for optimization variable)
%              - problem.ub (upper bound for optimization variable)
%
% Outputs:
%    q - semi-definite cost matrix
%    f - linear cost vector
%    a, blc, buc - two-sided constraint blc <= a*x <= buc
%    blx, bux - two-sided constraint blx <= x <= bux
%    param - empty []
%    cmd - minimization (enforced)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linprog2mosek, CORAquadprog

% Authors:       Mark Wetzlinger
% Written:       16-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only for quadprog, but we don't check
%   strcmp(problem.solver,'quadprog')
% for speed reasons

% note: all variables need to be here (can be empty)!
q = problem.H;
if ~isempty(problem.f)
    c = reshape(problem.f,[],1);
else
    c = zeros(length(q),0);
end
a = [problem.Aineq; problem.Aeq];
blc = [-Inf(size(problem.Aineq,1),1); problem.beq];
buc = [problem.bineq; problem.beq];
blx = problem.lb;
bux = problem.ub;

% mskqpopt expects a to be at least zeros(0,n)...
if isempty(a)
    a = sparse(0,length(q));
end

% no param
param = [];

% minimization, no display
cmd = 'minimize echo(0)';

% ------------------------------ END OF CODE ------------------------------
