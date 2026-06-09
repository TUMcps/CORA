function [val, x] = supportFunc_linprog(dir, Aineq, bineq, Aeq, beq, lb, ub, type)
% supportFunc_linprog - computes the support function via linear programming
%
% Syntax:
%    [val,x] = supportFunc_linprog(dir,Aineq,bineq,Aeq,beq,lb,ub,type)
%
% Inputs:
%    dir - cost function direction (column vector); the support function
%          value equals dir'*x for the upper bound
%    Aineq - inequality constraint matrix (Aineq * x <= bineq)
%    bineq - inequality constraint offset
%    Aeq - equality constraint matrix (Aeq * x == beq)
%    beq - equality constraint offset
%    lb - lower bound for optimization variable
%    ub - upper bound for optimization variable
%    type - 'upper', 'lower', or 'range'
%
% Outputs:
%    val - support function value (scalar for 'upper'/'lower', interval for
%          'range')
%    x - optimizer (column vector for 'upper'/'lower', [x_lower x_upper]
%        matrix for 'range')
%
% Other m-files required: CORAlinprog
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, polytope/supportFunc_, conZonotope/supportFunc_

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       26-March-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set up base problem struct (constraints stay the same for all types)
problem.Aineq = Aineq;
problem.bineq = bineq;
problem.Aeq = Aeq;
problem.beq = beq;
problem.lb = lb;
problem.ub = ub;

if strcmp(type,'lower')
    [val, x] = aux_solve(problem, dir);

elseif strcmp(type,'upper')
    [val, x] = aux_solve(problem, -dir);
    val = -val;

elseif strcmp(type,'range')
    % solve lower bound first
    [val_lower, x_lower, exitflag] = aux_solve(problem, dir);

    if exitflag == -2
        % primal infeasible -> empty set
        val = interval.empty(1);
        x = [];
        return
    end

    % use lower bound solution as initial point for upper bound
    problem.x0 = x_lower;
    [val_upper, x_upper] = aux_solve(problem, -dir);
    val_upper = -val_upper;

    % combine results
    val = interval(val_lower, val_upper);
    x = [x_lower, x_upper];
end

end


% Auxiliary functions -----------------------------------------------------

function [val, x, exitflag] = aux_solve(problem, f)
% solve LP: minimize f'*x subject to constraints

problem.f = f';

[x, val, exitflag] = CORAlinprog(problem);

if exitflag == -2
    % primal infeasible -> empty set
    val = Inf; x = [];
elseif exitflag == -3
    % unbounded
    val = -Inf;
    x = -sign(f) .* Inf(length(f),1);
elseif exitflag ~= 1
    throw(CORAerror('CORA:solverIssue'));
end

end

% ------------------------------ END OF CODE ------------------------------
