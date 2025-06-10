function empty = priv_representsa_emptySet(A,b,Ae,be,n,tol)
% priv_representsa_emptySet - determines if a polytope is empty
%
% Description:
%       A*x <= b
%    is infeasible if
%        min  b'*y
%        s.t. A'*y =  0,
%                y >= 0
%    is feasible. If problem is either unbounded (below since minimization)
%    or b'*y < 0, the polytope is empty. If the problem is infeasible or
%    b'*y == 0, the polytope is not empty.
%
% Syntax:
%    empty = priv_representsa_emptySet(A,b,Ae,be,n,tol)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    n - dimension of the polytope
%    tol - tolerance
%
% Outputs:
%    empty - true/false whether polytope is empty
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special case: 1D
if n == 1    
    % compute vertices (fast for 1D) and check whether empty
    V = priv_vertices_1D(A,b,Ae,be);
    empty = isempty(V);
    return
end

% quick check: no constraints
if isempty(b) && isempty(be)
    empty = false;
    return
end

% quick check: are there constraints of the form 0*x <= b with b < 0?
zero_rows = find(all(withinTol(A,0,tol),2));
if any(zero_rows) && any(b(zero_rows) < 0)
    empty = true;
    return
end


% number of inequality and equality constraints
nrConIneq = numel(b);
nrConEq = numel(be);

% normalize constraints
[A,b,Ae,be] = priv_normalizeConstraints(A,b,Ae,be,'A');

if ~isempty(Ae)
    % polytope can be already be declared empty if there is not point
    % that satisfies the set of equality constraints
    if rank(Ae) == min([nrConEq,n])
        % if the rank of the equality constraint matrix is smaller than
        % the minimum of the number of equality constraints and the
        % dimension, then there are still degrees of freedom
        
        % check for aligned vectors: pre-compute dot product
        dotprod_norm = Ae * Ae';
        
        for i=1:nrConEq
            % if there are aligned constraints, then they have to have
            % the same value in be (all but one of them are redundant),
            % otherwise the polytope is the empty set because there is
            % no point that can satisfy to parallel Ae*x = be at the
            % same time

            % check for aligned vectors
            alignedConstraints = all(...
                withinTol(Ae(i,:) - dotprod_norm(:,i)*Ae(i,:),...
                zeros(1,n),tol), 2);
        
            if nnz(alignedConstraints) > 1
                % at least two constraints are aligned
                if ~all(withinTol(be(alignedConstraints),be(i),tol))
                    % polytope is the empty set
                    empty = true;
                    return
                end
            end
        end

    end

    % rewrite equality constraints as inequality constraints
    [A,b] = priv_equalityToInequality(A,b,Ae,be);
    % update number of inequality constraints
    nrConIneq = length(b);
end


% solve the dual problem using linear programming

% set options for MATLAB linprog
persistent options
if isempty(options)
    options = optimoptions('linprog','display','off', ...
                       'OptimalityTolerance',1e-10);
end

problem.f = b';
problem.Aineq = [-eye(nrConIneq);A';-A'];
problem.bineq = zeros(nrConIneq+2*n,1);
problem.Aeq = [];
problem.beq = [];
problem.lb = [];
problem.ub = [];
problem.options = options;

% solve linear program (all inequalities)
[x,~,exitflag] = CORAlinprog(problem);

% read out solution
if exitflag == -2 || (exitflag > 0 && b'*x >= -1e-10)
    % if this problem is infeasible or if optimal objective value
    % is 0, the polytope is not empty
    empty = false;
elseif exitflag == -3 || (exitflag > 0 && b'*x < 1e-10)
    % if problem is unbounded (below since minimization) or
    % objective value is smaller zero, polytope is empty
    empty = true;
elseif exitflag == 0
    % max iterations exceeded; can not guarantee emptiness
    empty = false;
else
    throw(CORAerror('CORA:solverIssue','linprog'));
end

% ------------------------------ END OF CODE ------------------------------
