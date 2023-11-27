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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if property is set
if ~isempty(P.bounded.val)
    res = P.bounded.val;
    return
end

% if polytope has V-rep then it is bounded (1D might still have Inf)
if ~isempty(P.V.val) && ~any(any(isinf(P.V.val)))
    res = true;
    P.bounded.val = true;
    return
end

% dimension of ambient space
n = dim(P);

% 1D case
if n == 1
    % compute vertices (and save)
    V = vertices_(P,'lcon2vert');
    P.V.val = V;
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
    P.bounded.val = res;
    return
end

% quick check: if any column in the inequality and equality constraints is
% zero everywhere, this dimension is unbounded
if ~all(any([P.A; P.Ae],1))
    res = false;
    P.bounded.val = res;
    return
end

if ~isempty(P.fullDim.val) && P.fullDim.val
    % faster method if polytope is known to be non-degenerate
    res = aux_nondegenerate(P,n);
else
    res = aux_supportFunction(P,n);
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_nondegenerate(P,n)
% chedk via duality (only for degenerate polytopes!)

% the origin has to be contained in the (primal) polytope
if ~contains(P,zeros(n,1))
    % compute Chebyshev center and shift polytope
    c = center(P);
    if any(isnan(c))
        % set is unbounded (however, not every unbounded set yields a NaN
        % Chebyshev center!)
        res = false;
        P.bounded.val = false;
        P.emptySet.val = false;
        return
    elseif isempty(c)
        % set is empty
        res = true;
        P.bounded.val = true;
        P.emptySet.val = true;
        P.V.val = zeros(n,0);
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
P_ = polytope([P_.A; P_.Ae; -P_.Ae],[P_.b; P_.be; -P_.be]);

% normalize constraints to obtain dual polytope
%   A x <= 1  -->  P^* = conv(A^T)
P_ = normalizeConstraints(P_,'b');
% number of vertices
h = size(P_.A,1);

% check if the origin is contained in the dual polytope by linear
% programming:
%   min -t
%   s.t. A^T x = 0
%        for all i in {1,...,h}: t <= x_i
%        sum_i x_i = 1

% linear program options (only if MATLAB linprog is used)
persistent options
if isempty(options)
    options = optimoptions('linprog','display','off');
end

% variable vector: t in R, x in R^n
f = [-1; zeros(h,1)];

% inequality constraints: t <= x_i -> t - x_i <= 0
Aineq = [ones(h,1) -eye(h)];
bineq = zeros(h,1);

% equality constraints: A^T x = 0, sum_i x_i = 1
Aeq = [zeros(n,1), P_.A'; 0 ones(1,h)];
beq = [zeros(n,1); 1];

% solve linear program
[x,fval,exitflag] = linprog(f,Aineq,bineq,Aeq,beq,[],[],options);

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

% save to properties
P.bounded.val = res;

end

function res = aux_supportFunction(P,n)

    % loop over all dimensions
    for i=1:n
        for e=[-1,1]
            % check every axis-aligned direction
            dir = [zeros(i-1,1);e;zeros(n-i,1)];
            % evaluate support function
            val = supportFunc_(P,dir,'upper');
            if val == Inf
                % set is unbounded
                res = false;
                P.bounded.val = false;
                P.emptySet.val = false;
                return
            elseif val == -Inf
                % set is empty
                res = true;
                P.bounded.val = true;
                P.emptySet.val = true;
                return
            end
        end
    end

    % code reaches this part: set is neither unbounded nor empty
    res = true;
    P.bounded.val = true;
    P.emptySet.val = false;

end

% ------------------------------ END OF CODE ------------------------------
