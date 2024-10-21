function res = isequal(P,S,varargin)
% isequal - checks if a polytope is equal to another set or point
%
% Syntax:
%    res = isequal(P,S)
%    res = isequal(P,S,tol)
%
% Inputs:
%    P - polytope object 
%    S - contSet object 
%    tol - (optional) tolerance
%
% Outputs:
%    res - result of comparison
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([0 -1; 0 1;-1 0; 1 0; -1 -1],[2;3;2;3;2]);
%    res = isequal(P1,P2);
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

% Authors:       Mark Wetzlinger
% Written:       27-April-2023 (MW, migrated code from eq)
% Last update:   27-July-2023 (MW, handle 1D case)
%                16-December-2023 (MW, comparison to numeric)
%                02-December-2023 (MW, fix fully empty polytopes)
%                07-June-2024 (MW, quicker check for non-degenerate case)
%                14-July-2024 (MW, integrate comparison of V-polytopes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% set default values
tol = setDefaultValues(1e-6,varargin);

% check input arguments
inputArgsCheck({{P,'att',{'polytope','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% ensure that numeric is second input argument
[P,S] = reorderNumeric(P,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < P.precedence
    res = isequal(S,P,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(P,S,true)
    res = false;
    return
end

% handle fully empty polytope objects
if representsa_(P,'fullspace',0)
    res = representsa_(S,'fullspace',0);
    return
elseif representsa_(S,'fullspace',0)
    % if P were fullspace, we would have entered to if-branch above
    res = false;
    return
end

% general one-dimensional case
if dim(P) == 1
    res = aux_isequal_1D(P,S,tol); return
end

if isa(S,'polytope')
    % quick check: check emptiness, boundedness, and degeneracy
    if aux_diffProperties(P,S)
        res = false; return
    end

    % fast check if both are V-polytopes
    if P.isVRep.val && S.isVRep.val
        res = aux_isequal_Vpoly_Vpoly(P,S);
        return
    end

    % faster method if both are non-degenerate (only 2LPs to check)
    P_isnondeg = isFullDim(P);
    S_isnondeg = isFullDim(S);

    if P_isnondeg ~= S_isnondeg
        res = false;
    elseif P_isnondeg % S_isnondeg == true, too
        % both are non-degenerate
        res = aux_isequal_nondeg(P,S,tol);
    else
        % for degenerate polytopes, we check mutual containment
        res = contains_(P,S,'exact',tol) && contains_(S,P,'exact',tol);
    end
    return
end

% numeric
if isnumeric(S) && size(S,2) == 1
    % polytope must represent a single point
    if ~representsa_(P,'point',tol)
        res = false;
    else
        % compute that single point and check for equality up to tolerance
        V = vertices(P);
        res = all(withinTol(V,S,tol));
    end
    return
end

% check mutual containment for all other set representations
if isa(S,'contSet')
    res = contains_(P,S,'exact',tol) && contains_(S,P,'exact',tol);
    return
end

throw(CORAerror('CORA:noops',P,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_diffProperties(P1,P2)
% returns true if two polytopes P1 and P2 have the different properties (if
% given), e.g., P1 is bounded but P2 is unbounded; this immediately ensures
% that the polytopes cannot be equal!
% if properties unknown, we return false

% assume same properties (or unknown)
res = false;

% boundedness
if ~isempty(P1.bounded.val) && ~isempty(P2.bounded.val) ...
        && xor(P1.bounded.val,P2.bounded.val)
    res = true; return
end

% emptiness
if ~isempty(P1.emptySet.val) && ~isempty(P2.emptySet.val) ...
        && xor(P1.emptySet.val,P2.emptySet.val)
    res = true; return
end

% degeneracy
if ~isempty(P1.fullDim.val) && ~isempty(P2.fullDim.val) ...
        && xor(P1.fullDim.val,P2.fullDim.val)
    res = true; return
end

end

function res = aux_isequal_1D(P,S,tol)
% special method for 1D polytopes using vertex enumeration and supporting
% unbounded V-polytopes

% compute vertices (fast)
V1 = vertices_(P,'lcon2vert');
if isnumeric(S)
    V2 = S;
else
    V2 = vertices(S);
end

% compare
if xor(any(isinf(V1)),any(isinf(V2)))
    res = false;
elseif size(V1,2) ~= size(V2,2) 
    res = false;
elseif any(isinf(V1))
    res = all(V1 == V2);
else
    res = compareMatrices(V1,V2,tol);
end

end

function res = aux_isequal_Vpoly_Vpoly(P,S)
% check for set equality of two V-polytopes using [1, (7)]

% probably fastest to compute the minimal representation using the
% Quickhull algorithm instead of checking for the mutual containment of all
% points using LPs
tol = 1e-10;
P_V = priv_compact_V(P.V_.val,tol);
S_V = priv_compact_V(S.V_.val,tol);

% the matrices storing the vertices must be equal up to permutation
res = compareMatrices(P_V, S_V, tol, "equal");

end

function res = aux_isequal_nondeg(P,S,tol)
% quicker check for set equality between two non-degenerate polytopes

% ensure that both have H representation
[P_A,P_b,P_Ae,P_be] = constraints(P);
[S_A,S_b,S_Ae,S_be] = constraints(S);

% row-wise normalization
[P_A,P_b,P_Ae,P_be] = priv_normalizeConstraints(P_A,P_b,P_Ae,P_be,'A');
[S_A,S_b,S_Ae,S_be] = priv_normalizeConstraints(S_A,S_b,S_Ae,S_be,'A');

% rewrite equalities as inequalities (otherwise we'd have to cross check)
% and unify with offset vector
[P_A,P_b] = priv_equalityToInequality(P_A,P_b,P_Ae,P_be);
[S_A,S_b] = priv_equalityToInequality(S_A,S_b,S_Ae,S_be);
P_Ab = [P_A,P_b];
S_Ab = [S_A,S_b];

% loop over all halfspaces of P
% save indices of matched rows in S_Ab
idxMatched = false(size(S_Ab,1),1);
nrCon_P = size(P_Ab,1);
for j=1:nrCon_P
    % check if j-th row is also in other matrix
    idxInS = all(withinTol(S_Ab - P_Ab(j,:), 0, tol),2);

    % if there is no match with any halfspace in S, then the halfspace must
    % be redundant, otherwise the polytopes are different
    if ~any( idxInS ) && aux_isNonredundantConstraint(P_Ab,j,tol)
        res = false; return;
    end

    % update list of matched halfspaces in S
    idxMatched = idxMatched | idxInS;
end

% ...repeat same process as above switching P and S
nrCon_S = size(S_Ab,1);
for j=1:nrCon_S
    % loop only of those halfspaces in S that have not yet been matched
    if ~idxMatched(j)
        % check if j-th row is also in other matrix
        idxInP = all(withinTol(S_Ab(j,:) - P_Ab, 0, tol),2);
    
        % if there is no match with any halfspace in P, then the halfspace
        % must be redundant, otherwise the polytopes are different
        if ~any( idxInP ) && aux_isNonredundantConstraint(S_Ab,j,tol)
            res = false; return;
        end
    end
end

% if the code reaches here, then every halfspace in either polytope either
% matches a halfspace of the other polytope or it is redundant, hence the
% polytopes are equal
res = true;

end

function res = aux_isNonredundantConstraint(P_Ab, j, tol)
% compute support function in the direction of the j-th normal vector of
% the polytope where the j-th constraint has been removed
% ...since this would require instantiating many polytopes, we copy the
% code from the supportFunc_ for now
dir = P_Ab(j,1:end-1)';

problem.f = -dir';
problem.Aineq = [P_Ab(1:j-1,1:end-1);P_Ab(j+1:end,1:end-1)];
problem.bineq = [P_Ab(1:j-1,end); P_Ab(j+1:end,end)];
problem.Aeq = [];
problem.beq = [];
problem.lb = [];
problem.ub = [];

% solve linear program
[~,val,exitflag] = CORAlinprog(problem);
val = -val;

if exitflag == -3
    % unbounded
    val = Inf;
elseif exitflag == -2
    % infeasible -> empty set
    val = -Inf;
elseif exitflag ~= 1
    throw(CORAerror('CORA:solverIssue'));
end

% compare the value of the support function with the offset from the
% constraint: if the value is smaller or equal up to tolerance, then the
% j-th constraint is determined to be redundant
% (note: can also handle Inf values!)
res = val > P_Ab(j,end) && ~withinTol(val, P_Ab(j,end), tol);

end

% ------------------------------ END OF CODE ------------------------------
