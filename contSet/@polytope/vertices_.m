function V = vertices_(P,method,varargin)
% vertices_ - computes the vertices of a polytope; this function also
%    serves as a getter function for the property 'V' of a polytope object
%
% Syntax:
%    V = vertices_(P,method)
%
% Inputs:
%    P - polytope object
%    method - 'lcon2vert' (default, duality method),
%             'comb' (all combinations of n inequalities)
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

% if polytope has V-representation, return it
if P.isVRep.val
    V = P.V;
    return
end

% dimension
n = dim(P);

% check if polytope is known to be empty
if ~isempty(P.emptySet.val) && P.emptySet.val
    V = zeros(n,0);
    P.isVRep.val = true;
    return
end

% 1D case quick
if n == 1
    V = aux_vertices_1D(P);
    P.V_.val = V;
    P.isVRep.val = true;
    return
end

% compute Chebyshev center to detect unbounded/empty cases (only one linear
% program, we will generally require several below anyway)
c = center(P);
if any(isnan(c))
    % polytope is unbounded
    P.emptySet.val = false;
    P.bounded.val = false;
    throw(CORAerror('CORA:notSupported',...
                    'Vertex enumeration requires a bounded polytope.'));
elseif isempty(c)
    % polytope is empty
    V = zeros(n,0);
    P.V_.val = V;
    P.isVRep.val = true;
    P.minVRep.val = true;
    P.emptySet.val = true;
    P.bounded.val = true;
    P.fullDim.val = false;
    return
end


% currently three different methods supported: we potentially fall back to
% 'lcon2vert' from 'cdd', hence the structure below
if strcmp(method,'cdd')
    [V,method] = aux_vertices_cdd(P,c);
end

if strcmp(method,'lcon2vert')
    V = aux_vertices_lcon2vert(P,n,c);
elseif strcmp(method,'comb')
    % note: this method computes a superset containing the vertices;
    % caution: the convex hull of these points is an outer approximation of
    % the original polytope! (implemented for debugging purposes)
    V = aux_vertices_comb(P);
end

% set properties
P.V_.val = V;
P.isVRep.val = true;
P.minVRep.val = true; %...?
% emptiness has been check above
P.emptySet.val = false;
% unbounded cases would have already thrown an error
P.bounded.val = true;
% determine degeneracy via SVD
if size(V,2) <= size(V,1)
    % full-dimensionality requires at least n+1 vertices
    P.fullDim.val = false;
else
    [~,S,~] = svd(V - mean(V,2));
    P.fullDim.val = n == nnz(~withinTol(S,0,1e-12));
end

end


% Auxiliary functions -----------------------------------------------------

function V = aux_vertices_1D(P)
% simplified function for 1D polytopes    

    % compute minimal representation
    P = compact_(P,'all',1e-9);
    P = normalizeConstraints(P,'A');

    % if there is no point, P is already empty
    if P.emptySet.val
        V = []; return
    end

    if ~isempty(P.A_.val)

        % check boundedness from below
        Aisminus1 = withinTol(P.A_.val,-1);
        if any(Aisminus1)
            % bounded from below
            V = -P.b_.val(Aisminus1);
        else
            % unbounded toward -Inf
            V = -Inf;
        end
    
        % check boundedness from above
        Ais1 = withinTol(P.A_.val,1);
        if any(Ais1)
            % bounded from above (add only if not a duplicate)
            if ~withinTol(V,P.b_.val(Ais1))
                V = [V, P.b_.val(Ais1)];
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
        V = P.be_.val / P.Ae_.val;
        % equality representation -> not full-dimensional
        P.fullDim.val = false;
        % emptiness would have been detected above
        P.emptySet.val = false;

    else
        throw(CORAerror('CORA:specialError',...
            'Error in vertex computation of 1D polytope.'));

    end

end

function [V,method] = aux_vertices_cdd(P,c)
% vertex enumeration using cdd (third-party)

% cdd requires shift by Chebyshev center according to mpt toolbox
P_centered = P - c;

% call cddmex... returns struct with vertices as a field
try
    s = cddmex('extreme', ...
        struct('A',[P_centered.Ae_.val;P_centered.A_.val],...
               'B',[P_centered.be_.val;P_centered.b_.val],...
               'lin',1:size(P_centered.Ae_.val,1)));
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

function V = aux_vertices_lcon2vert(P,n,c)
% vertex enumeration using 'lcon2vert' method

% lcon2vert only supports non-degenerate polytopes, otherwise throws an
% error: in that case, we compute the basis of the affine hull and compute
% the vertices in that subspace
try
    % vertex enumeration algorithm
    v = lcon2vert(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,[],0);
    % transpose so that each column is a vertex
    V = v';

catch ME
    % check for subspaces
    [~,S] = isFullDim(P);

    % is polytope just a single point?
    if isempty(S)
        V = c;
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
        P_subspace = polytope(P.A_.val*S,P.b_.val-P.A_.val*c,...
            P.Ae_.val*S,P.be_.val-P.Ae_.val*c);
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
% an error and we would have exited in the catch-branch above)
% -> check unboundedness via duality (quick since we have V representation)
if ~isBounded(P)
    P.bounded.val = false;
    throw(CORAerror('CORA:notSupported',...
            'Vertex computation requires a bounded polytope.'));
end

end

function V = aux_vertices_comb(P)
% simple vertex enumeration algorithm: this function returns a set of
% vertices that contains the true vertices; however, the minimal vertices
% are not the convex hull of the computed vertices

% check if polytope is unbounded
if ~isBounded(P)
    throw(CORAerror('CORA:notSupported',...
        'Vertex computation requires a bounded polytope.'));
end

if ~isempty(P.be_.val)
    % rewrite as inequality constraints
    A = [P.A_.val; P.Ae_.val; -P.Ae_.val];
    b = [P.b_.val; P.be_.val; -P.be_.val];
    P = polytope(A,b);
end

% normalize rows
P = normalizeConstraints(P,'A');

% minimal H-representation
P = compact_(P,'all',1e-9);

% read out constraints
A = P.A_.val;
b = P.b_.val;

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

% ------------------------------ END OF CODE ------------------------------
