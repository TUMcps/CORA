function [x,fval,exitflag,output,lambda] = CORAquadprog(problem)
% CORAquadprog - evaluates a quadratic program for MATLAB and MOSEK syntax;
%    we need this function for the following reason: MOSEK does not support
%    the syntax
%       quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options)
%    as it cannot deal with the options struct, if given. Thus, we use
%       quadprog(problem)
%    for the MATLAB call and reformulate the variables for MOSEK here.
%
% Syntax:
%    [x,fval,exitflag,output,lambda] = CORAquadprog(problem)
%
% Inputs:
%    problem - quadratic program definition, with mandatory fields
%              - problem.H (quadratic cost function min f*x)
%              - problem.f (linear cost function min f*x)
%              - problem.A (inequality constraint Aineq * x <= bineq)
%              - problem.b (inequality constraint Aineq * x <= bineq)
%              - problem.Aeq (equality constraint Aeq * x == beq)
%              - problem.beq (equality constraint Aeq * x == beq)
%              - problem.lb (lower bound for optimization variable)
%              - problem.ub (upper bound for optimization variable)
%
% Outputs:
%    x - minimizer
%    fval - minimal objective value
%    exitflag - status of linear program
%    output - further output of MATLAB linprog
%    lambda - further output of MATLAB linprog
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadprog2mosek

% Authors:       Mark Wetzlinger
% Written:       18-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if MOSEK is installed
persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end

if isMosek
    % convert problem struct to MOSEK syntax
    [q,c,a,blc,buc,blx,bux,param,cmd] = quadprog2mosek(problem);

    % call MOSEK
    res = mskqpopt(q,c,a,blc,buc,blx,bux,param,cmd);

    % convert to MATLAB linprog outputs
    if strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
        % feasible
        exitflag = 1;

        x = [];
        if isfield(res,'sol')
            if isfield(res.sol,'itr')
                x = res.sol.itr.xx;
            else
                x = res.sol.bas.xx;
            end
        end
        % objective value: evaluate quadratic cost function
        fval = 0.5*x'*q*x + c'*x;

    else
        % either infeasible, unbounded, or unknown error
        x = [];
        fval = [];
        
        if strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
            % infeasible
            exitflag = -2;
        elseif strcmp(res.sol.itr.prosta,'DUAL_INFEASIBLE')
            % unbounded
            exitflag = -3;
        else
            throw(CORAerror('CORA:solverIssue'));
        end
    end

    % set other output arguments
    if nargout >= 4
        output = [];
    end
    if nargout >= 5
        lambda = [];
    end

else
    % quadratic program options
    if ~isfield(problem,'options')
        persistent options
        if isempty(options)
            options = optimoptions('quadprog','display','off');
        end
        problem.options = options;
    end

    problem.solver = 'quadprog';

    % call MATLAB quadprog (with struct syntax)
    [x,fval,exitflag,output,lambda] = quadprog(problem);

end

% ------------------------------ END OF CODE ------------------------------
