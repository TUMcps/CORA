function res = isBounded(P)
% isBounded - determines whether a polytope is bounded
%
% Syntax:
%    res = isBounded(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    res - true/false whether the polytope is bounded
%
% Example: 
%    A = [2 1; 1 3; -1 2; -4 1];
%    b = ones(4,1);
%    P = polytope(A,b);
%    res = isBounded(P)
%
% Reference:
%    [1] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       31-May-2022
% Last update:   30-November-2022 (MW, add 1D case)
%                08-December-2022 (MW, quick exit in nD case)
%                27-July-2023 (MW, fix 1D case)
%                25-November-2023 (MW, faster method for degenerate polytopes)
%                04-June-2024 (MW, use simplex for support function evals)
% Last revision: 10-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% check if property is set
if ~isempty(P.bounded.val)
    res = P.bounded.val;
    return
end

% check whether V- or H-representation given
if P.isVRep.val
    % --- V representation
    res = aux_isBounded_Vpoly(P);

else
    % --- H representation

    % dimension of ambient space
    n = dim(P);
    
    % 1D case
    if n == 1
        res = aux_isBounded_1D_Hpoly(P);
    else
    
        % quick check: if any column in the inequality and equality
        % constraints is zero everywhere, this dimension is unbounded
        if ~all(any([P.A_.val; P.Ae_.val],1))
            res = false;
        elseif ~isempty(P.fullDim.val) && P.fullDim.val
            % faster method if polytope is known to be non-degenerate
            res = aux_isBounded_nD_Hpoly_nondegenerate(P,n);
        else
            % standard method
            res = aux_isBounded_nD_Hpoly(P,n);
        end

    end

end

% save the set property (only done once, namely, here!)
P.bounded.val = res;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isBounded_Vpoly(P)

% if polytope has V-rep then it is bounded, unless it is in 1D, then
% we still have to check for Inf
res = dim(P) > 1 || ~any(any(isinf(P.V_.val)));

end

function res = aux_isBounded_1D_Hpoly(P)

% compute vertices
V = vertices_(P,'lcon2vert');

% empty vertices or non-Inf vertices -> bounded
if isempty(V)
    res = true;
    P.emptySet.val = true;
elseif any(isinf(V))
    res = false;
    P.emptySet.val = false;
else
    res = true;
    P.emptySet.val = false;
end

end

function res = aux_isBounded_nD_Hpoly_nondegenerate(P,n)
% check via duality (only for degenerate polytopes!)

% the origin has to be contained in the (primal) polytope
if ~contains(P,zeros(n,1))
    % compute Chebyshev center and shift polytope
    c = center(P);
    if any(isnan(c))
        % set is unbounded (however, not every unbounded set yields a NaN
        % Chebyshev center!)
        res = false;
        P.emptySet.val = false;
        return
    elseif isempty(c)
        % set is empty
        res = true;
        P.emptySet.val = true;
        P.V = zeros(n,0);
        P.isVRep = true;
        P.minVRep.val = true;
        return
    end
    % set is not empty (otherwise c would have been empty)
    P.emptySet.val = false;

    % shift set by center
    P_ = P - c;
else
    % we need to reassign P due to the handling of the properties
    P_ = P;
end

% rewrite equality constraints as pairwise inequality constraints
[A,b] = priv_equalityToInequality(P_.A_.val,P_.b_.val,P_.Ae_.val,P_.be_.val);

% normalize constraints to obtain dual polytope
%   A x <= 1  -->  P^* = conv(A^T)
[A,b] = priv_normalizeConstraints(A,b,[],[],'b');
% number of vertices
h = size(A,1);

% check if the origin is contained in the dual polytope by linear
% programming:
%   min -t
%   s.t. A^T x = 0
%        for all i in {1,...,h}: t <= x_i
%        sum_i x_i = 1

% variable vector: t in R, x in R^n
problem.f = [-1; zeros(h,1)];

% inequality constraints: t <= x_i -> t - x_i <= 0
problem.Aineq = [ones(h,1) -eye(h)];
problem.bineq = zeros(h,1);

% equality constraints: A^T x = 0, sum_i x_i = 1
problem.Aeq = [zeros(n,1), A'; 0 ones(1,h)];
problem.beq = [zeros(n,1); 1];

problem.lb = [];
problem.ub = [];

% solve linear program
[x,fval,exitflag] = CORAlinprog(problem);

if exitflag >= 0
    % solution could be computed
    if fval > 0 || withinTol(fval,0,1e-8)
        % origin is not contained in dual polytope -> primal polytope is
        % unbounded
        res = false;
    else
        res = true;
    end
else
    % TODO: handle other cases...
    throw(CORAerror('CORA:solverIssue'));
end

end

function res = aux_isBounded_nD_Hpoly(P,n)
% check if all support functions in the directions of an n-dimensional
% simplex are bounded; directions are [I_n, -1_n]

    % loop over columns of identity matrix
    for i=1:n
        % check every axis-aligned direction
        direction = unitvector(i,n);
        % evaluate support function
        val = supportFunc_(P,direction,'upper');
        if val == Inf
            % set is unbounded
            res = false;
            P.emptySet.val = false;
            return
        elseif val == -Inf
            % set is empty
            res = true;
            P.emptySet.val = true;
            return
        end
    end
    % check -1_n
    direction = -ones(n,1);
    % evaluate support function
    val = supportFunc_(P,direction,'upper');
    if val == Inf
        % set is unbounded
        res = false;
        P.emptySet.val = false;
        return
    elseif val == -Inf
        % set is empty
        res = true;
        P.emptySet.val = true;
        return
    end

    % code reaches this part: set is neither unbounded nor empty
    res = true;
    P.emptySet.val = false;

end

% ------------------------------ END OF CODE ------------------------------
