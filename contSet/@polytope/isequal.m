function res = isequal(P,S,varargin)
% isequal - check ifs two polytopes are equal
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default values
tol = setDefaultValues(1e-6,varargin);

% check dimensions
try
    equalDimCheck(P,S);
catch ME
    % P and S have a different ambient dimension
    res = false; return
end

% handle fully empty polytope objects
if representsa_(P,'fullspace',0)
    res = representsa_(S,'fullspace',0); return;
elseif representsa_(S,'fullspace',0)
    % if P were fullspace, we would have entered to if-branch above
    res = false; return
end

% one-dimensional case
if dim(P) == 1
    res = aux_1D(P,S,tol); return
end

% quick check for polytopes: check emptiness, boundedness, and degeneracy
if isa(S,'polytope') && aux_diffProperties(P,S)
    res = false; return
end

if isa(S,'polytope')
    % check degeneracy: faster method if both are non-degenerate
    P_isnondeg = isFullDim(P);
    S_isnondeg = isFullDim(S);

    if P_isnondeg ~= S_isnondeg
        res = false; return;
    elseif P_isnondeg
        % both are non-degenerate
        res = aux_isequal_nondeg(P,S,tol);
    else
        % for degenerate polytope, we use mutual containment
        res = contains_(P,S,'exact',tol) && contains_(S,P,'exact',tol);
    end

elseif isnumeric(S)
    % avoid comparison with matrices
    if size(S,2) > 1
        throw(CORAerror('CORA:noops',P,S));
    end

    % S is a single point
    if ~representsa_(P,'point',tol)
        res = false;
    else
        % compute vertex
        V = vertices(P);
        res = all(withinTol(V,S,tol));
    end

else
    % check mutual containment for all other set representations
    res = contains_(P,S,'exact',tol) && contains_(S,P,'exact',tol);

end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_nondeg(P,S,tol)
% quicker check for set equality between two non-degenerate polytopes

% row-wise normalization
P = normalizeConstraints(P,'A');
S = normalizeConstraints(S,'A');

% rewrite equalities as inequalities (otherwise we'd have to cross check)
% and unify with offset vector
P_Ab = [[P.A; P.Ae; -P.Ae], [P.b; P.be; P.be]];
S_Ab = [[S.A; S.Ae; -S.Ae], [S.b; S.be; S.be]];

% loop over all halfspaces of P
% save indices of matched rows in S_Ab
idxMatched = false(size(S_Ab,1),1);
nrCon_P = size(P_Ab,1);
for j=1:nrCon_P
    % check if j-th row is also in other matrix
    idxInS = all(withinTol(S_Ab - P_Ab(j,:), 0, tol),2);

    % if there is no match with any halfspace in S, then the halfspace must
    % be redundant, otherwise the polytopes are different
    if ~any( idxInS ) && aux_isNonredundant(P_Ab,j,tol)
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
        if ~any( idxInP ) && aux_isNonredundant(S_Ab,j,tol)
            res = false; return;
        end
    end
end

% if the code reaches here, then every halfspace in either polytope either
% matches a halfspace of the other polytope or it is redundant, hence the
% polytopes are equal
res = true;

end

function res = aux_isNonredundant(P_Ab, j, tol)
% compute support function in the direction of the j-th normal vector of
% the polytope where the j-th constraint has been removed
% ...since this would require instantiating many polytopes, we copy the
% code from the supportFunc_ for now
dir = P_Ab(j,1:end-1)';

persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end

persistent options
if isempty(options)
    options = optimoptions('linprog','display','off');
end

% solve linear program
if isMosek

    % number of constraints
    nrCon = size(P_Ab,1);

    % rewrite for MOSEK syntax
    c = dir';
    a = [P_Ab(1:j-1,1:end-1); P_Ab(j+1:end,1:end-1)];
    blc = -Inf(nrCon-1,1);
    buc = [P_Ab(1:j-1,end); P_Ab(j+1:end,end)];
    
    % call MOSEK
    res = msklpopt(c,a,blc,buc,[],[],[],'maximize echo(0)');

    % read out value of support function and support vector
    if strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
        val = res.sol.itr.dobjval;
    else
        % either unbounded or infeasible
        if strcmp(res.sol.itr.prosta,'DUAL_INFEASIBLE')
            % unbounded
            val = Inf;
        elseif strcmp(res.sol.itr.prosta,'PRIMAL_INFEASIBLE')
            % infeasible -> empty set
            val = -Inf;
        end
    end

else

    problem.f = -dir';
    problem.Aineq = [P_Ab(1:j-1,1:end-1);P_Ab(j+1:end,1:end-1)];
    problem.bineq = [P_Ab(1:j-1,end); P_Ab(j+1:end,end)];
    problem.solver = 'linprog';
    problem.options = options;

    % solve linear program
    [~,val,exitflag] = linprog(problem);
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
end

% compare the value of the support function with the offset from the
% constraint: if the value is smaller or equal up to tolerance, then the
% j-th constraint is determined to be redundant
% (note: can also handle Inf values!)
res = val > P_Ab(j,end) && ~withinTol(val, P_Ab(j,end), tol);

end

function res = aux_1D(P,S,tol)

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

function res = aux_diffProperties(P1,P2)
% returns true if two polytopes P1 and P2 have the different properties (if
% given), e.g., P1 is bounded but P2 is unbounded
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

% ------------------------------ END OF CODE ------------------------------
