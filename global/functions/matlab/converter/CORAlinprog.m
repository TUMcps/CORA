function [x,fval,exitflag,output,lambda] = CORAlinprog(problem)
% CORAlinprog - evaluates a linear program for MATLAB and MOSEK syntax;
%    we need this function for the following reasons: MATLAB's linprog
%    function could previously be called by
%       linprog(f,A,b,Ae,be,lb,ub,x0,options)
%    which was recently changed to
%       linprog(f,A,b,Ae,be,lb,ub,options)
%    forcing us to use the alternative syntax
%       linprog(problem)
%    where problem is a struct. However, the MOSEK overload of linprog
%    cannot handle that syntax, so we use this wrapper, instead.
%
% Syntax:
%    [x,fval,exitflag,output,lambda] = CORAlinprog(problem)
%
% Inputs:
%    problem - linear program definition, with mandatory fields
%              - problem.f (cost function min f*x)
%              - problem.Aineq (inequality constraint Aineq * x <= bineq)
%              - problem.bineq (inequality constraint Aineq * x <= bineq)
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-July-2024
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
    [f,a,blc,buc,blx,bux,param,cmd] = linprog2mosek(problem);

    % call MOSEK
    res = msklpopt(f,a,blc,buc,blx,bux,param,cmd);

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
        % objective value: evaluate cost function
        fval = f'*x;

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
    % linear program options
    if ~isfield(problem,'options')
        persistent options
        if isempty(options)
            options = optimoptions('linprog','display','off');
        end
        problem.options = options;
    end

    problem.solver = 'linprog';

    % call MATLAB linprog (with struct syntax)
    [x,fval,exitflag,output,lambda] = linprog(problem);

end

% ------------------------------ END OF CODE ------------------------------
