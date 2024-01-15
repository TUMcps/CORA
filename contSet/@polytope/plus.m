function P_out = plus(P,S)
% plus - overloaded '+' operator for the Minkowski addition of two
%    polytopes or a polytope with a vector
%
% Syntax:
%    P_out = plus(P,S)
%
% Inputs:
%    P - polytope object or numerical vector
%    S - polytope object or numerical vector
%
% Outputs:
%    P_out - polytope after Minkowski addition
%
% Example:
%    A = [2 1; -1 1; -2 -3; 0 -4; 2 -1]; b = ones(5,1);
%    P = polytope(A,b);
%    v = [2;1];
%    res = P + v;
%
%    figure; hold on;
%    plot(P);
%    plot(res,[1,2],'r');
%
% References: MPT-toolbox https://www.mpt3.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Viktor Kotsev
% Written:       20-June-2022
% Last update:   25-October-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% sort arguments (Minkowski addition is commutative)
[P,S] = findClassArg(P,S,'polytope');

% check dimensions
equalDimCheck(P,S);

% dimension
n = dim(P);

% polytope + polytope
if isa(S,"polytope")
    %check for empty polytopes
    if representsa(P, 'emptySet')
        P_out = S; return
    elseif representsa(S, 'emptySet')
        P_out = P; return
    end
    % check if one of the sets is the origin
    if representsa(P,'origin')
        P_out = S; return
    elseif representsa(S,'origin')
        P_out = P; return
    end

    % both polytopes have vertex representation
    if ~isempty(P.V.val) && ~isempty(S.V.val)
        % read out vertices
        V1 = P.V.val; V2 = S.V.val;
        numV1 = size(V1,2); numV2 = size(V2,2);

        % init resulting vertices
        V = zeros(n,numV1*numV2);

        % add each vertex to each vertex
        for i = 1:numV2
			V(:,((i-1)*numV1+1):i*numV1,:) = bsxfun(@plus, V1, V2(:,i));
        end
        
        % init resulting polytope
        P_out = polytope(V);
        return
    end

    % all-zero matrices
    PZ = zeros(size(P.A,1),size(S.A,2));
    PZe = zeros(size(P.Ae,1),size(S.Ae,2));

    % lift inequalities and equalities to higher-dimensional space
    A = [S.A -S.A; PZ P.A];
    b = [S.b; P.b];
    Ae = [S.Ae -S.Ae; PZe P.Ae];
    be = [S.be; P.be];

    % project resulting polytope onto original dimensions
    P_highdim = polytope(A,b,Ae,be);
    P_out = project(P_highdim,1:n);

elseif isnumeric(S)
    % avoid addition with matrices
    if size(S,2) > 1
        throw(CORAerror('CORA:noops',P,S));
    end

    % shift offsets
    b = P.b + P.A*S;
    be = P.be + P.Ae*S;

    % init output polytope
    P_out = polytope(P.A,b,P.Ae,be);

    % shift vertices if V representation is given
    if ~isempty(P.V.val)
        P_out.V.val = P.V.val + S;
    end

    % assign properties (same as P)
    P_out = copyProperties(P,P_out,'noV');

elseif isa(S,'zonotope') || isa(S,'interval') || ...
    isa(S,'conZonotope') || isa(S,'zonoBundle')
    S = polytope(S);
    P_out = P + S;

elseif isa(S,'polyZonotope')
    P_out = S + P;

else
    throw(CORAerror('CORA:noops',P,S));
    
end

% set properties

% If both polytopes are bounded, then sum is also bounded
if isa(P_out,"polytope") && (~isempty(P.bounded.val) && P.bounded.val) ...
    && isa(S,"polytope") && (~isempty(S.bounded.val) && S.bounded.val)
    P_out.bounded.val = true;
end

% If one of the polytopes is unbounded, then sum is also unbounded
if isa(P_out,"polytope") && (~isempty(P.bounded.val) && ~P.bounded.val) ...
    ||  (isa(S,"polytope") && ~isempty(S.bounded.val) && ~S.bounded.val)
    P_out.bounded.val = false;
end

% If one of the polytopes is fully dimensional, then sum is also fully dimensional
if isa(P_out,"polytope") && (~isempty(P.fullDim.val) && P.fullDim.val) ...
    ||  (isa(S,"polytope") && ~isempty(S.fullDim.val) && S.fullDim.val)
    P_out.fullDim.val = true;
end

end

% ------------------------------ END OF CODE ------------------------------
