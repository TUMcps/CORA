function V = vertices_(P,method,varargin)
% vertices_ - computes the vertices of a polytope
%
% Syntax:
%    V = vertices_(P,method)
%
% Inputs:
%    P - polytope object
%    method - 'lcon2vert' (default, duality method),
%             'comb' (all combinations of n inequalities, then removal of
%                     redundancies)
%             'cdd' (double descriptor method, requires cddmex)
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
%                21-November-2023 (MW, improve comb algorithm)
%                24-November-2023 (MW, implement cdd method)
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

% compute Chebyshev center to detect unbounded/empty cases
c = center(P);
if any(isnan(c))
    % polytope is unbounded
    P.emptySet.val = false;
    P.bounded.val = false;
    throw(CORAerror('CORA:notSupported',...
                    'Vertex computation requires a bounded polytope.'));
elseif isempty(c)
    % polytope is empty
    V = zeros(n,0);
    P.V.val = V;
    P.minVRep.val = true;
    P.emptySet.val = true;
    P.bounded.val = true;
    P.fullDim.val = false;
    return
end

if strcmp(method,'cdd')
    % cdd requires shift by Chebyshev center according to mpt toolbox
    P_centered = P - c;

    % call cddmex... returns struct with vertices as a field
    try
        s = cddmex('extreme', ...
            struct('A',[P_centered.Ae;P_centered.A],...
                   'B',[P_centered.be;P_centered.b],...
                   'lin',1:size(P_centered.Ae,1)));
        % if s.R is non-empty, there are vertex rays -> unbounded
        if ~isempty(s.R)
            P_centered.emptySet.val = false;
            P_centered.bounded.val = false;
            throw(CORAerror('CORA:notSupported',...
                    'Vertex computation requires a bounded polytope.'));
        end
        V = s.V' + c;
    catch ME
        if strcmp(ME.identifier,'CORA:notSupported')
            rethrow(ME);
        end
        % fallback option: switch to lcon2vert method below...
        method = 'lcon2vert';
    end
end

if strcmp(method,'lcon2vert')
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
            % within the original polytope P (we use the origin since the
            % center has been subtracted above)
            c = center(P);
            P_subspace = polytope(P.A*S,P.b-P.A*c,P.Ae*S,P.be-P.Ae*c);
            % compute vertices for y in subspace
            V = vertices_(P_subspace,'lcon2vert');
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

    % polytope is not degenerate (otherwise lcon2vert would have thrown
    % an error and we would have exited above)
    % -> check unboundedness via duality (quick)
    P.fullDim.val = true;
    if ~isBounded(P)
        P.bounded.val = false;
        throw(CORAerror('CORA:notSupported',...
                'Vertex computation requires a bounded polytope.'));
    end
end

if strcmp(method,'comb')
    % check if polytope is unbounded
    if ~isBounded(P)
        throw(CORAerror('CORA:notSupported',...
            'Vertex computation requires a bounded polytope.'));
    end

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
P.minVRep.val = true; %...?
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
% some combinations don't produce a vertex or the vertex is outside the 
% polytope... use logical indexing at the end to avoid resizing the matrix
% that contains all vertices
idxKeep = true(1,nrComb);

% toggle warning, since some intersection points will be -/+Inf
warOrig = warning;
warning('off','all');

% loop over all combinations
for i=1:nrComb
    % take n halfspaces from A
    A_ = A(comb(i,:),:);
    % if full-rank, then we have exactly one intersection point
    if rank(A_,1e-8) < n
        idxKeep(i) = false;
        continue
    end

    % compute intersection point of n halfspaces taken from A
    V(:,i) = A_ \ b(comb(i,:));

    % check if vertex is contained in polytope
    val = A*V(:,i);
    if ~all( val < b | withinTol(val,b,1e-8) )
        idxKeep(i) = false;
        continue
    end
    
    % check if vertex is a duplicate
    if i > 1
        dist = vecnorm(V(:,1:i-1) - V(:,i));
        idxKeep(i) = ~any(withinTol(dist(idxKeep(1:i-1)),0,1e-14));
    end

end

warning(warOrig);

% remove vertices at indices where there was no computation, or the
% computed vertex is outside of the polytope
V = V(:,idxKeep);

% remove vertices with Inf/NaN values
V = V(:,all(isfinite(V),1));

end

function V = aux_1D(P)
% simplified function for 1D polytopes    

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
