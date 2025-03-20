function c = center(P,varargin)
% center - computes the Chebyshev center of a polytope
%    note: polytope in vertex representation are converted to halfspace
%    representation, see [1], which is potentially time-consuming
%    Use method 'avg' for average of vertices
%
% Syntax:
%    c = center(P)
%
% Inputs:
%    P - polytope object
%    method - 'chebyshev', 'avg' (for v-polytope)
%
% Outputs:
%    c - Chebyshev center of the polytope
%
% Example:
%    P = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    c = center(P);
%    
% Reference:
%    [1] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger, Tobias Ladner
% Written:       28-March-2022
% Last update:   14-December-2022 (MW, add call to MOSEK)
%                27-July-2023 (MW, add fast method for 1D)
%                02-January-2024 (MW, fix fully empty polytopes)
%                31-October-2024 (TL, added methods, added 'avg')
%                25-February-2025 (TL, unconstrained polytope returns origin)
% Last revision: 12-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% parse input
method = setDefaultValues({'chebyshev'},varargin);
allowedMethods = {'chebyshev','avg'};
inputArgsCheck({ ...
    {P,'att','polytope'}; ...
    {method,'str',allowedMethods}; ...
})

% read out dimension
n = dim(P);

% fullspace/empty case
if representsa_(P,'fullspace',0)
    % return origin; consistent with fullspace/center
    c = zeros(n,1); return;
elseif ~isempty(P.emptySet.val) && P.emptySet.val
    % return empty
    c = double.empty(n,0); return;
end

% fast and simple computation for 1D
if n == 1
    c = aux_center_1D(P); return
end

% switch method
switch method
    case 'chebyshev'
        c = aux_center_chebyshev(P);
    case 'avg'
        c = aux_center_avg(P);
    otherwise
        throw(CORAerror('CORA:wrongValue','second',allowedMethods))
end

end


% Auxiliary functions -----------------------------------------------------

function c = aux_center_1D(P)
% special method for 1D polytopes

% compute vertices
V = vertices_(P,'lcon2vert');
if isempty(V)
    % empty set
    c = zeros(1,0);
elseif size(V,2) == 1
    % only one vertex, which is also the center
    c = V;
elseif any(isinf(V))
    % unbounded
    c = NaN;
else
    % bounded -> two vertices, take mean
    c = mean(V);
end

end

function c = aux_center_chebyshev(P)

% read out dimension
n = dim(P);

% check whether there are only equalities: allows to avoid the linear
% program from below (faster)
if isempty(P.A_.val) && ~isempty(P.Ae_.val)
    c = aux_center_only_equalityConstraints(P,n);
    return
end

% general method: compute Chebyshev center via linear program; to this end,
% we require the halfspace representation
constraints(P);
c = aux_center_LP(P,n);

end

function c = aux_center_avg(P)
    % compute vertices
    V = vertices(P);

    % compute mean
    c = mean(V,2);
end

function c = aux_center_only_equalityConstraints(P,n)
% three outcomes: unbounded, single point, infeasible

% minimal halfspace representation: if two constraints are aligned and
% cannot be fulfilled at the same time, an empty polytope is returned
[~,~,Ae,be] = priv_normalizeConstraints([],[],P.Ae_.val,P.be_.val,'A');
[Ae,be,empty] = priv_compact_alignedEq(Ae,be,1e-12);

% check if emptyness has been determined during the computation of the
% minimal representation
if empty
    c = double.empty(n,0);
    return
end

% all constraints now are linearly independent, hence the relation of 
% system dimension and number of constraints determines the solution
if size(Ae,1) < n
    % underdetermined -> unbounded
    c = NaN(n,1);
elseif size(Ae,1) > n
    % overdetermined -> no solution
    c = double.empty(n,0);
else
    % same number of constraints as system dimension -> single point
    c = Ae \ be;
end

end

function c = aux_center_LP(P,n)
% linear program for the computation of the Chebyshev center

% dimension and number of (in)equalities
nrEq = size(P.Ae_.val,1);

% 2-Norm of each row
A_norm = sqrt(sum(P.A_.val.^2,2));

% extend inequality and equality constraints by one column
A_ext = [P.A_.val, A_norm];
Ae_ext = [P.Ae_.val, zeros(nrEq,1)];

% cost function for linear program: minimize 2-norm of constraints
f = [zeros(n,1); -1];

% init linprog struct
problem.f = f;
problem.Aineq = A_ext;
problem.bineq = P.b_.val;
problem.Aeq = Ae_ext;
problem.beq = P.be_.val;
problem.lb = [-Inf(n,1);0];
problem.ub = [];

% solve LP
[c,val,exitflag] = CORAlinprog(problem);

if exitflag == 1
    % truncate solution
    c = c(1:n);
elseif exitflag == -2
    % set is empty
    c = double.empty(n,0);
elseif exitflag == -3
    % unbounded
    c = NaN(n,1);
elseif exitflag < 0
    throw(CORAerror('CORA:solverIssue'));
end

end

% ------------------------------ END OF CODE ------------------------------
