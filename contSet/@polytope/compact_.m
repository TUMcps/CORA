function P_out = compact_(P,type,tol,varargin)
% compact_ - Computes the irredundant H/V-representation of a polytope;
%    inequality constraints: the main idea is to compare the offset of a
%    given constraint with the support function of the polytope without
%    that constraints; specialized methods for 1D and 2D polytopes
%
% Syntax:
%    P_out = compact_(P)
%    P_out = compact_(P,type)
%    P_out = compact_(P,type,tol)
%
% Inputs:
%    P - polytope object
%    type - (optional) set of constraints that should be reduced
%           'all' (default): inequality and equality constraints
%           'A': only inequality constraints
%           'Ae': only equality constraints
%           'aligned': check only for constraints with same orientation
%           'zeros': check for 0*x <= a and 0*x = a constraints
%           'V': remove redundancies in vertex representation
%           'AtoAe': rewrite pairwise inequalities as equalities
%
% Outputs:
%    P_out - polytope object in reduced or minimal representation
%
% Example:
%    A = [1 1; 1 0; 1 1; -1 0; 0 -1; 1 1; 1 2];
%    b = [7; 3; 5; -1; -1; 6; 20];
%    P = polytope(A,b);
%    P_ = compact(P);
%
%    figure; hold on;
%    plot(P);
%    plot(P_,[1,2],'r--');
%
% Reference: MPT-Toolbox https://www.mpt3.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/compact

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       23-May-2022
% Last update:   08-December-2022 (MW, use logical indexing to eliminate inequalities in sequential tests)
%                09-December-2022 (MW, pre-processing using orthogonality and alignment criteria)
%                12-December-2022 (MW, fix faulty alignment criterion)
%                14-December-2022 (MW, add support for MOSEK)
%                04-April-2023 (MW, use persistent variables)
%                15-January-2024 (TL, bug fix single infeasible equality constraint)
%                25-September-2024 (MW, add variant for inequalities to equalities)
% Last revision: 28-May-2023 (MW, increase readability via aux_ functions, integrate code from removeRedundancies)
%                31-July-2023 (MW, rename 'compact_', integrate deleteZeros, integrate minVRep)

% ------------------------------ BEGIN CODE -------------------------------

% vertex representation
if strcmp(type,'V')
    if ~P.isVRep.val
        throw(CORAerror('CORA:notSupported',...
            'Vertex representation not provided.'));
    end
    V = priv_compact_V(compact_(P,'lcon2vert').V,tol);
    P_out = polytope(P);
    P_out.V_.val = V;
    P_out.minVRep.val = true;
    return
end

%%% what to do with below?
% quick exits:
% if length(P_out.b_.val) == 1 && isempty(P_out.be_.val)
%     % only one inequality, no equality given -> already minimal
%     P_out = aux_oneHalfspace(P_out,tol);
%     P = aux_copyProperties(P,P_out);
%     return
% elseif isempty(P_out.b_.val) && length(P_out.be_.val) == 1
%     % only one equality, no inequality given -> already minimal
%     P_out.minHRep.val = true;
%     P = aux_copyProperties(P,P_out);
%     return
% end
% % quick check for two inequality constraints
% if length(P_out.b_.val) == 2 && isempty(P_out.Ae_.val)
%     P_out = aux_twoHalfspaces(P_out);
%     P_out.minHRep.val = true;
%     % check whether object has become fully empty
%     P = aux_copyProperties(P,P_out);
%     return
% end


% normalized halfspace representation
n = dim(P);
[A,b,Ae,be] = constraints(P);
[A,b,Ae,be] = priv_normalizeConstraints(A,b,Ae,be,'A');
% unless we use the most general method, we do not obtain the minimal
% halfspace representation
minHRep = false;
empty = false;

switch type
    case 'zeros'
        [A,b,Ae,be,empty] = priv_compact_zeros(A,b,Ae,be,tol);

    case 'aligned'
        % remove aligned halfspaces
        [A,b] = priv_compact_alignedIneq(A,b,tol);
        [Ae,be,empty] = priv_compact_alignedEq(Ae,be,tol);

    case 'AtoAe'
        % convert inequality to equality constraints
        [A,b,Ae,be] = priv_compact_toEquality(A,b,Ae,be,tol);

    case 'A'
        [A,b,Ae,be] = priv_compact_toEquality(A,b,Ae,be,tol);
        [A,b] = priv_compact_alignedIneq(A,b,tol);
        [A,b,~,~,empty] = priv_compact_nD(A,b,[],[],n,tol);

    case 'Ae'
        [Ae,be,empty] = priv_compact_alignedEq(Ae,be,tol);

    case 'all'
        % put the logic into a separate function since we cannot 'break'
        % the case statement, but return prematurely from a function
        [A,b,Ae,be,empty,minHRep] = priv_compact_all(A,b,Ae,be,n,tol);
end

% emptiness detected
if empty
    P_out = polytope.empty(n);
    P = priv_copyProperties(P_out,P,'all');
    return
end

% if no constraints are left, we have an unbounded set (fullspace)
if isempty(A) && isempty(Ae)
    P_out = polytope.Inf(n);
    P = priv_copyProperties(P_out,P,'all');
    return
end

% instantiate polytope and copy properties
P_out = polytope(A,b,Ae,be);
P_out = priv_copyProperties(P,P_out,'all');
if minHRep
    P_out.minHRep.val = true;
end

% ------------------------------ END OF CODE ------------------------------
