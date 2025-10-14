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
%    Since the dual-simplex algorithm sometimes returns exitflag = -9 for
%    problems which do have a solution, we have an automatic switch to the
%    interior-point solver in that case. Note that we use the dual-simplex
%    algorithm by default since it has shown to be more accurate.
%
% Syntax:
%    [x,fval,exitflag,output,lambda] = CORAlinprog(problem)
%
% Inputs:
%    problem - linear program definition, with fields
%              - problem.f (cost function min f*x)
%              - problem.Aineq (inequality constraint Aineq * x <= bineq)
%              - problem.bineq (inequality constraint Aineq * x <= bineq)
%              - problem.Aeq (equality constraint Aeq * x == beq)
%              - problem.beq (equality constraint Aeq * x == beq)
%              - problem.lb (lower bound for optimization variable)
%              - problem.ub (upper bound for optimization variable)
%              - problem.x0 (initial point)
%              where all numeric values should be of type double.
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
% Last update:   04-October-2024 (MW, switch to interior-point in MATLAB call)
%                09-October-2024 (TL, compatibility >=R2024a)
%                29-October-2024 (TL, problem fields should be doubles)
%                30-May-2026 (MP, Included new option with higher precision)
%                27-August-2025 (TL, compatibility >=R2025a)
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
        % Sometimes, the interior-point solver does a great job, but
        % doesn't know about it, and outputs an 'UNKNOWN', then calling
        % upon the might of the simplex solver to get to an optimal
        % solution. Let us check whether this is possible:
        try % We are not certain that the simplex algorithm has indeed been
            % used.
            prosta = res.sol.bas.prosta;
        catch ME
            prosta = '';
        end

        % Let us check what the simplex solver yields:
        if strcmp(res.sol.bas.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
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
    % (either use user-defined one or automatically determine best one)
    if isfield(problem,'options')
        options = {problem.options};
    else
        options = {};
        
        % default algorithm: 'dual-simplex'
        persistent optionsDualSimplexDefault
        if isempty(optionsDualSimplexDefault)
            optionsDualSimplexDefault = optimoptions('linprog','Algorithm','dual-simplex','display','off');
        end
        options = [options {optionsDualSimplexDefault}];

        if ~isMATLABReleaseOlderThan('R2024a') && isMATLABReleaseOlderThan('R2025a')
            % set of options with lower tolerances in case of very small
            % values occuring in the LP
            % Contained within this 'if' block since it leads to problems
            % with old matlab versions
            optionsHighTol = optimoptions('linprog','Algorithm','dual-simplex','display','off','ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
            options = [options {optionsHighTol}];
            
            % matlab updated the 'dual-simplex' algorithm in R2024a ('dual-simplex-highs').
            % However, the old algorithm provides solutions in some situations when the new one does not (exitflag=0)..
            % We add the old algorithm here as fallback option.
            if isMATLABReleaseOlderThan('R2025a') % only available until <R2025a...
                persistent optionsDualSimplexLegacy
                if isempty(optionsDualSimplexLegacy)
                    w = warning; warning off; % don't show legacy warning
                    optionsDualSimplexLegacy = optimoptions('linprog','Algorithm','dual-simplex-legacy','display','off');
                    warning(w);
                end
                options = [options {optionsDualSimplexLegacy}];
            end
        end
    end

    % test all options
    for i=1:numel(options)

        problem.options = options{i};
        algorithm = problem.options.Algorithm;
        problem.solver = 'linprog';
    
        % call MATLAB linprog (with struct syntax)
        [x,fval,exitflag,output,lambda] = aux_MATLABlinprog(problem,algorithm);
    
        % in some cases, the dual-simplex algorithm loses feasibility (leading
        % to exitflag = -9), so we switch to the interior-point solver and try
        % again
        if exitflag == -9 && startsWith(algorithm,'dual-simplex')
            algorithm = 'interior-point';
            [x,fval,exitflag,output,lambda] = aux_MATLABlinprog(problem,algorithm);
        end

        if  ~isMATLABReleaseOlderThan('R2025a') && exitflag == -2 && startsWith(algorithm,'dual-simplex')
            % TL: this sometimes happens since R2025a, I contacted the support
            % and they told me they are working on it... in the meantime, we
            % should try the interior point method ...
            algorithm = 'interior-point';
            vargoutIP = cell(1,5);
            [vargoutIP{:}] = aux_MATLABlinprog(problem,algorithm);
            if vargoutIP{3} > 0
                % overwrite results if a solution is found
                [x,fval,exitflag,output,lambda] = vargoutIP{:};
            else
                % also try legacy
                algorithm = 'interior-point-legacy';
                vargoutIPL = cell(1,5);
                [vargoutIPL{:}] = aux_MATLABlinprog(problem,algorithm);
                if vargoutIPL{3} > 0
                    % overwrite results if a solution is found
                    [x,fval,exitflag,output,lambda] = vargoutIPL{:};
                end
            end
        end

        if exitflag > 0
            % some solution is found, don't test other options
            break;
        end

        if exitflag == -3
            % problem is unbounded, we assume there are no false-positives
            % here and exit
            break
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [x,fval,exitflag,output,lambda] = aux_MATLABlinprog(problem,algorithm)

% call MATLAB linprog
w = warning; warning off; % don't show legacy warning
try
    problem.options.Algorithm = algorithm; % set the solver
    [x,fval,exitflag,output,lambda] = linprog(problem);

catch ME
    if strcmp(ME.identifier,'optim:linprog:NonDoubleInput')
       % convert relevant fields of problem to double
       % (done here for efficiency)

       % show warning
       warning on;
       CORAwarning("CORA:solver",'Not all given fields for the linear program are of type double. Re-trying with converted values.')
       warning off;

       % convert fields
       fields = {'f','Aineq','bineq','Aeq','beq','lb','ub','x0'};
       for i=1:numel(fields)
           field = fields{i};
           % if field exists
           if isfield(problem,field) && isnumeric(problem.(field))
               % convert field to double
               problem.(field) = double(problem.(field));
           end
       end

       % re-call linprog
       [x,fval,exitflag,output,lambda] = linprog(problem);
    elseif strcmp(ME.identifier,'optim:linprog:UnknownError') && startsWith(algorithm,'dual-simplex')
        % TL: this sometimes happens since R2025a, I contacted the support
        % and they told me they are working on it... in the meantime, we
        % should try the interior point method.
        try
            algorithm = 'interior-point';
            [x,fval,exitflag,output,lambda] = aux_MATLABlinprog(problem,algorithm);
        catch
           rethrow(ME) 
        end
    else
        % rethrow all other exceptions
        rethrow(ME)
    end
end
% reset warning level
warning(w);

end

% ------------------------------ END OF CODE ------------------------------
