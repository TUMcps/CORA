function [res,subspace] = isFullDim(P)
% isFullDim - checks if the dimension of the affine hull of a polytope is
%    equal to the dimension of its ambient space; additionally, one can
%    obtain a basis of the subspace in which the polytope is contained
%
% Syntax:
%    res = isFullDim(P)
%    [res,subspace] = isFullDim(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    res - true/false
%    subspace - (optional) Returns a set of orthogonal unit vectors
%               x_1,...,x_k such that P is strictly contained in
%               center(P)+span(x_1,...,x_k)
%               (here, 'strictly' means that k is minimal).
%               Note that if P is just a point, subspace=[].
%
% Example:
%    P1 = polytope([1 0;0 1;-1 0; 0 -1],[1;1;0;0]);
%    P2 = polytope([1 0;0 1;-1 0; 0 -1],[1;0;0;0]);
%    isFullDim(P1)
%    isFullDim(P2)
%
% Reference:
%    [1] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Authors:       Niklas Kochdumper, Viktor Kotsev, Mark Wetzlinger, Adrian Kulmburg
% Written:       02-January-2020 
% Last update:   08-August-2022
%                05-December-2022 (VK, new method)
%                14-December-2022 (MW, add call to MOSEK)
%                25-May-2023 (AK, added method to compute subspace)
%                28-May-2023 (MW, add quick check for pairwise annihilation)
%                27-July-2023 (MW, add 1D method)
%                03-May-2024 (MW, fix subspace computation for V-rep)
% Last revision: 10-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% check if fullDim property is already set
if ~isempty(P.fullDim.val)
    res = P.fullDim.val;
    if nargout == 1
        return
    else
        % compute subspace below...
    end
end

% check whether V- or H-representation given
if P.isVRep.val
    % --- V representation
    [res, subspace] = aux_isFullDim_Vpoly(P,nargout);

else
    % --- H representation
    [res, subspace] = aux_isFullDim_Hpoly(P,nargout);

end

% save the set property (only done once, namely, here!)
P.fullDim.val = res;

end


% Auxiliary functions -----------------------------------------------------

function [res,subspace] = aux_isFullDim_Vpoly(P,num_argout)
% compute degeneracy acc. to [1, (17)] and subspace acc. to [1, (20)]

% set tolerance, read out dimension
tol = 1e-10;
n = dim(P);

% compute rank shifted by mean
V_shifted = P.V_.val - mean(P.V_.val,2);
rankV = rank(V_shifted,tol);
% compare to ambient dimension
res = rankV == n;

% compute basis of affine hull
if num_argout == 2
    if res
        subspace = eye(n);
    else
        [Q,R] = qr(V_shifted);
        subspace = Q(:,1:rankV);
    end
else
    % dummy value for subspace (not used further)
    subspace = NaN;
end

end

function [res,subspace] = aux_isFullDim_Hpoly(P,num_argout)

if dim(P) == 1
    % 1D case
    [res,subspace] = aux_isFullDim_1D_Hpoly(P);
else
    % >=2D case
    if num_argout < 2
        % If the user only wants to know whether the polytope is
        % degenerate, we can check this rapidly using linear
        % programming, see [1, (18)]
        res = aux_isFullDim_nD_Hpoly_nosubspace(P);
        % put a dummy value here (will not be used further)
        subspace = NaN;
    
    elseif num_argout == 2
        % If on the other hand, the user is interested in the subspace,
        % we have to do the computations from above in a more specific
        % way, see [1, Alg. 2]
        [res,subspace] = aux_isFullDim_nD_H_subspace(P);
        
    end
end

end

function [res,subspace] = aux_isFullDim_1D_Hpoly(P)
% special method for 1D polytopes (vertices fast to compute)

% compute vertices
V = vertices_(P,'lcon2vert');

% go over cases...
if isempty(V)
    % empty set
    res = false; subspace = [];
    P.emptySet.val = true;
    P.bounded.val = true;
elseif size(V,2) == 1
    % single vertex
    res = false; subspace = [];
    P.emptySet.val = false;
    P.bounded.val = true;
elseif any(isinf(V))
    % unbounded
    res = true; subspace = 1;
    P.emptySet.val = false;
    P.bounded.val = false;
else
    % bounded
    res = true; subspace = 1;
    P.emptySet.val = false;
    P.bounded.val = true;
end

end

function res = aux_isFullDim_nD_Hpoly_nosubspace(P)

% Existence of equality constraints implies degenerate case (set may
% also be empty)
if ~isempty(P.Ae_.val)
    res = false;
    return
end

% number of inequality constraints and dimension
[nrCon,n] = size(P.A_.val);

% normalize constraints
P_norm = normalizeConstraints(P,'A');

% quick check: any pair of inequalities of type
%   ax <= b, ax >= b_ (-ax <= -b_), with b_ >= b?
dotprod = tril(P_norm.A * P_norm.A');
idxAntiParallel = find(withinTol(dotprod,-1));
for i=1:length(idxAntiParallel)
    % read out matrix indexing from linear indexing
    [r,c] = ind2sub([nrCon,nrCon],idxAntiParallel(i));
    % check offset values, two perspectives:
    %  ax <= b1, -ax <= b2 <=> ax >= -b2  ->  -b2 >= b1?
    % -ax <= b1 <=> ax >= -b1, ax <= b2   ->  -b1 >= b2 <=> -b2 >= b1?
    % ...the one with the reversed sign has to be larger (or equal)!
    if -P_norm.b(c) > P_norm.b(r) || withinTol(-P_norm.b(c),P_norm.b(r))
        res = false;
        return
    end
end

% 2-Norm of each row (always 1 since we normalized the constraints)
A_norm = ones(nrCon,1);

% extend inequality and equality constraints by one column
A = [P.A_.val,A_norm];

% cost function for linear program: minimize 2-norm of constraints
f = [zeros(n,1); -1];

% init linprog struct
problem.f = f;
problem.Aineq = A;
problem.bineq = P.b_.val;
problem.Aeq = [];
problem.beq = [];
problem.lb = [-Inf(n,1);0];
problem.ub = [];

% Solve Linear Program
[~,r,exitflag] = CORAlinprog(problem);

if exitflag == 1
    % check if point is strictly interior
    if withinTol(0,r,1e-9)
        % r is 0 -> no ball fits in polytope -> degenerate
        res = false;
    else
        % -r < 0 => radius of ball greater than 0 that fits in
        % polytope -> non-degenerate
        res = true;
    end
elseif exitflag == -2
    % set is empty
    res = false;
    P.emptySet.val = true;
elseif exitflag == -3
    % numerical issues seem to come up only in the case of
    % full-dimensional unbounded polytopes
    res = true;
elseif exitflag < 0
    throw(CORAerror('CORA:solverIssue'));
end

end

function [res,subspace] = aux_isFullDim_nD_H_subspace(P)
% since we have to evaluate many linear programs anyway, we can quickly
% check if the polytope is completely empty -> also degenerate
if representsa_(P,'emptySet',1e-10)
    res = false;
    subspace = [];
    return
end

% read out dimension
n = dim(P);

% Let's first check whether the polytope is empty, by searching for one
% point that is contained in P; we don't call center(P) since the function
% below also supports unbounded sets in a more useful way (i.e., not
% returning NaN), but it fails if the polytope is only the origin, which is
% why we check whether the origin is contained beforehand
if contains_(P,zeros(n,1),'exact',1e-5)
    P_iter = polytope(P);
else
    x0 = aux_maxNormPerpendicularPolytope(P,zeros([n,0]));
    if isempty(x0)
        res = false;
        P.emptySet.val = true;
        subspace = [];
        return
    end
    
    % Now, assuming x0 exists, we translate P by -x0 so that we can
    % guarantee that the resulting polytope contains the origin:
    P_iter = P - x0;
end

% Setup the list of vectors we seek
subspace = zeros([n,0]); % Need to have this form to avoid bugs

for i = 1:n
    % Search for vectors in P_iter that are perpendicular to the ones
    % we have found so far
    x_iter = aux_maxNormPerpendicularPolytope(P_iter,subspace);
    % If the solution is x_iter=0 (within same tolerance as used for
    % solving LPs), then this means that there are no non-trivial
    % vectors that are perpendicular to those we have found, and so we
    % can stop the iteration
    if norm(x_iter,Inf) <= 1e-8
        break;
    end
    
    x_iter_unit = x_iter / norm(x_iter);
    
    % Now, it suffices to add x_iter to the list
    subspace = [subspace,x_iter_unit];
    
end

% Now that we have constructed the subspace we need, it is time to
% check what dimension it has; only if it is full dimensional, is the
% polytope non-degenerate
k = size(subspace,2);
if k == n
    res = true;
    % If that is the case, we can assume that subspace is the entire
    % space, and so a much easier ONB would be the following
    subspace = eye(n);
else
    res = false;
end

end

function x = aux_maxNormPerpendicularPolytope(P,X)
% For a list of vectors X = [x_1,...,x_k], solve the linear program
%   max ||x||_oo
%   s.t. x + Xw \in P,
%        forall i: x_i'*x = 0
%
% This is equivalent to
%   max_y max_x y'*x,
%   s.t. x + Xw \in P,
%        forall i: x_i'*x = 0,
%        and where y is iterated over all basis vectors +-e_i.
    
% read out dimension of polytope
n = dim(P);
% store maximum objective value and maximizer
maximum = 0;
maximizer = [];

% loop over all dimensions (for y)
for i=1:n
    % compute maximum for y = +e_i
    y = zeros(n,1);
    y(i) = 1;
    [res, x] = aux_firstMaximum(P,y,X);

    % save maximizer if objective value is larger
    if -res > maximum
        maximum = -res;
        maximizer = x;
    end
    
    % compute maximum for y = -e_i
    y = zeros(n,1);
    y(i) = -1;
    [res, x] = aux_firstMaximum(P,y,X);

    % save maximizer if objective value is larger
    if -res > maximum
        maximum = -res;
        maximizer = x;
    end
end

% return maximizer
x = maximizer;

end

function [res, x] = aux_firstMaximum(P,y,X)
    
% Extend to have x + Xw \in P
P_A_extended = [P.A_.val, P.A_.val*X];
P_Aeq_extended = [P.Ae_.val, P.Ae_.val*X];
y_extended = [y; zeros([size(X,2),1])];

% add constraints forall i: x_i'*x=0
P_Aeq_extended = [P_Aeq_extended; X' zeros([size(X,2),size(X,2)])];
P_beq_extended = [P.be_.val; zeros(size(X,2),1)];

% init linprog struct
problem.f = -y_extended;
problem.Aineq = P_A_extended;
problem.bineq = P.b_.val;
problem.Aeq = P_Aeq_extended;
problem.beq = P_beq_extended;
problem.lb = [];
problem.ub = [];

[x,res,exitflag] = CORAlinprog(problem);

% If the problem is unbounded, we need to add a constraint, e.g.,
%   y'*x = 1
% in order to find a good direction for x.
if exitflag == -3
    problem.Aeq = [P_Aeq_extended; y_extended'];
    problem.beq = [P_beq_extended; 1];
    [x,~] = CORAlinprog(problem);
    % set the objective value manually to -Inf to force updating the
    % maximizer (see calling function)
    res = -Inf;
end

x = x(1:dim(P));

end

% ------------------------------ END OF CODE ------------------------------
