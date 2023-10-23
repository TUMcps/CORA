function V = vertices_(P,method,varargin)
% vertices_ - computes the vertices of a polytope
%
% Syntax:
%    V = vertices_(P,method)
%
% Inputs:
%    P - polytope object
%    method - 'lcon2vert' (default),
%             'comb' (all combinations of n inequalities)
%
% Outputs:
%    V - vertices
%
% Example: 
%    A = [1 0 -1 0 1; 0 1 0 -1 1]';
%    b = [3; 2; 3; 2; 1];   
%    P = polytope(A,b);
%
%    V = vertices(P);
%
%    figure; hold on;
%    plot(P);
%    plot(V(1,:),V(2,:),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices, zonotope/vertices

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       01-February-2011
% Last update:   12-June-2015
%                12-August-2016
%                09-May-2018 (NK, changed algorithm for vertex computation)
%                29-June-2018 (MA, all polytopes considered)
%                13-December-2022 (MW, simple vertex enumeration algorithm)
%                30-May-2023 (MW, support degenerate cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% computation of minimal H-representation is contained in lcon2vert
% infeasibility check (isempty-like function) also contained in lcon2vert

% dimension
n = dim(P);

% check if polytope is known to be empty
if ~isempty(P.emptySet.val) && P.emptySet.val
    V = zeros(n,0); return
end

% if polytope has V-representation, return it
if ~isempty(P.V.val)
    if n == 1
        % reduce to minimal representation while we're at it
        temp = [min(P.V.val),max(P.V.val)];
        if ~any(isinf(temp)) && withinTol(temp(1),temp(2))
            temp = temp(1);
        end
        P.V.val = temp;
    end
    V = P.V.val; return
end

% 1D case quick
if n == 1
    V = aux_1D(P);
    P.V.val = V;
    return
end

if strcmp(method,'lcon2vert')

    % isempty is just one linear program (there are several in lcon2vert)
    if representsa(P,'emptySet')
        V = zeros(n,0);
        P.V.val = V;
        P.minVRep.val = true;
        P.emptySet.val = true;
        P.bounded.val = true;
        P.fullDim.val = false;
        return
    end
    
    try
        % vertex enumeration algorithm
        v = lcon2vert(P.A,P.b,P.Ae,P.be,[],0);
        % transpose so that each column is a vertex
        V = v';
    catch ME
        % check for subspaces
        [~,S] = isFullDim(P);

        % is polytope just a point?
        if isempty(S)
            V = center(P);
            P.V.val = V;
            P.minVRep.val = true;
            P.emptySet.val = false;
            P.fullDim.val = false; % no zero-dimensional sets
            P.bounded.val = true;
            return
        end

        % polytope is non-degenerate in subspace
        if size(S,2) < n
            % we plug in
            %    x = S*y + c in Ax <= b, resulting in A*(Sy) <= b - Ac,
            % where y is of the subspace dimension and c is any point
            % within the original polytope P (we use the center-function)
            c = center(P);
            P_ = polytope(P.A*S,P.b - P.A*c);
            % compute vertices for y in subspace
            V = vertices_(P_,'lcon2vert');
            % map back via x = S*y + c
            V = S*V + c;
            return
        elseif size(S,2) == n
            % check if polytope is unbounded
            if ~isBounded(P)
                throw(CORAerror('CORA:notSupported',...
                    'Vertex computation requires a bounded polytope.'));
            end
        end
        rethrow(ME);
    end

elseif strcmp(method,'comb')
    if ~isempty(P.be)
        % rewrite as inequality constraints
        A = [P.A; P.Ae; -P.Ae];
        b = [P.b; P.be; -P.be];
        P = polytope(A,b);
    end
    V = aux_simpleVertexEnum(P);
end

% set hidden properties
P.V.val = V;
P.minVRep = true; %...?
P.emptySet.val = false;
P.bounded.val = true;
% determine degeneracy via SVD
if size(V,2) <= size(V,1)
    % full-dimensionality requires at least n+1 vertices
    P.fullDim.val = false;
else
    [~,S,~] = svd(V);
    P.fullDim.val = size(V,1) == nnz(~withinTol(S,0,1e-12));
end

end


% Auxiliary functions -----------------------------------------------------

function V = aux_simpleVertexEnum(P)
% simple vertex enumeration algorithm: this function returns a set of
% vertices that contains the true vertices; however, the minimal vertices
% are not the convex hull of the computed vertices

% normalize rows
P = normalizeConstraints(P,'A');

% minimal H-representation
P = compact_(P,'all',1e-9);

% read out constraints
A = P.A;
b = P.b;

% number of constraints and dimension
[nrCon,n] = size(A);

% combinator does not work if n > nrCon, i.e., degenerate cases
if n > nrCon
    throw(CORAerror('CORA:notSupported',...
        'Method ''type'' does not support degenerate cases.'));
end

% all possible combinations of n constraints
nrComb = nchoosek(nrCon,n);

% throw error if computational effort too high
if nrComb > 10000
    throw(CORAerror('CORA:specialError','Too many combinations.'));
end

% list all combinations
comb = combinator(nrCon,n,'c');
% init vertices
V = zeros(n,nrComb);

% toggle warning, since some intersection points will be -/+Inf
warOrig = warning;
warning('off','all');

% loop over all combinations
for i=1:nrComb
    % compute intersection point of n halfspaces taken from A
    V(:,i) = A(comb(i,:),:) \ b(comb(i,:));
end

warning(warOrig);

% remove vertices with Inf/NaN values
V = V(:,all(isfinite(V),1));

end

function V = aux_1D(P)
    
    % compute minimal representation
    P = compact_(P,'all',1e-9);
    P = normalizeConstraints(P,'A');

    % if there is no point, P is already empty
    if P.emptySet.val
        V = []; return
    end

    if ~isempty(P.A)

        % check boundedness from below
        Aisminus1 = withinTol(P.A,-1);
        if any(Aisminus1)
            % bounded from below
            V = -P.b(Aisminus1);
        else
            % unbounded toward -Inf
            V = -Inf;
        end
    
        % check boundedness from above
        Ais1 = withinTol(P.A,1);
        if any(Ais1)
            % bounded from above (add only if not a duplicate)
            if ~withinTol(V,P.b(Ais1))
                V = [V, P.b(Ais1)];
            end
        else
            % unbounded toward +Inf
            V = [V, Inf];
        end

        % set boundedness, emptiness
        P.emptySet.val = false;
        P.bounded.val = any(isinf(V));
        P.fullDim.val = size(V,2) > 1;

    elseif ~isempty(P.Ae)
        % due to minHRep call above, we should only have one equality here
        V = P.be / P.Ae;
        % equality representation -> not full-dimensional
        P.fullDim.val = false;
        % emptiness would have been detected above
        P.emptySet.val = false;

    else
        throw(CORAerror('CORA:specialError',...
            'Error in vertex computation of 1D polytope.'));

    end

end

% ------------------------------ END OF CODE ------------------------------
