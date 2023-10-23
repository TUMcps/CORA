function [res,S] = representsa_(pZ,type,tol,varargin)
% representsa_ - checks if a polynomial zonotope can also be represented by
%    a different set, e.g., a special case
%
% Syntax:
%    res = representsa_(pZ,type,tol)
%    [res,S] = representsa_(pZ,type,tol)
%
% Inputs:
%    pZ - polyZonotope object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    S - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Mark Wetzlinger, Niklas Kochdumper, Victor Gassmann
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(pZ,type);
else
    [empty,res,S] = representsa_emptyObject(pZ,type);
end
if empty; return; end

% dimension
n = dim(pZ);

% init second output argument (covering all cases with res = false)
S = [];

% minimal representation
pZ = compact_(pZ,'all',eps);

% is polynomial zonotope a single point?
isPoint = isempty(pZ.G) && isempty(pZ.GI);

switch type
    case 'origin'
        res = ~isempty(pZ.c) && ( (isPoint && all(withinTol(pZ.c,0,tol))) ...
            || norm(interval(pZ)) <= tol);
        if nargout == 2 && res
            S = zeros(n,1);
        end

    case 'point'
        res = isPoint;
        if nargout == 2 && res
            S = pZ.c;
        end

    case 'capsule'
        % 1D, a point, or maximum one generator
        res = n == 1 || isPoint && size(pZ.G,2) + size(pZ.GI) <= 1;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from polyZonotope to ' type ' not supported.']));
        end

    case 'conHyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of polyZonotope to ' type ' not supported.']));
        
    case 'conPolyZono'
        res = true;
        if nargout == 2
            S = conPolyZono(pZ);
        end

    case 'conZonotope'
        res = n == 1 || aux_isZonotope(pZ,tol) || aux_isPolytope(pZ,tol);
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from polyZonotope to ' type ' not supported.']));
        end

    case 'ellipsoid'
        res = n == 1 || isPoint;
        if nargout == 2 && res
            S = ellipsoid(pZ);
        end

    case 'halfspace'
        % polynomial zonotopes are bounded
        res = false;

    case 'interval'
        res = n == 1 || (aux_isZonotope(pZ,tol) && representsa_(zonotope(pZ),'interval',tol));
        if nargout == 2 && res
            S = interval(pZ);
        end

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of polyZonotope to ' type ' not supported.']));

    case 'polytope'
        res = n == 1 || aux_isZonotope(pZ,tol) || aux_isPolytope(pZ,tol);
        if nargout == 2 && res
            S = polytope(pZ);
        end

    case 'polyZonotope'
        % obviously true
        res = true;
        if nargout == 2
            S = pZ;
        end

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        res = n == 1 || aux_isZonotope(pZ,tol) || aux_isPolytope(pZ,tol);
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from polyZonotope to ' type ' not supported.']));
        end

    case 'zonotope'
        res = n == 1 || aux_isZonotope(pZ,tol);
        if nargout == 2 && res
            S = zonotope(pZ);
        end

    case 'hyperplane'
        % polynomial zonotopes are bounded
        res = false;

    case 'parallelotope'
        res = n == 1 || (isempty(pZ.G) && representsa_(zonotope(pZ),'parallelotope',tol));
        if nargout == 2 && res
            S = zonotope(pZ);
        end

    case 'emptySet'
        % already handled in isemptyobject
        res = false;

    case 'fullspace'
        res = false;

end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isZonotope(pZ,tol)

    if isempty(pZ.G)
        res = true; return
    end

    res = false;
    
    % remove redundant exponent vectors
    [E,G] = removeRedundantExponents(pZ.E,pZ.G);
    
    % check matrix dimensions
    if size(E,1) ~= size(G,2)
        return;
    end
    
    % sort exponent matrix rows
    E = sortrows(E,'descend');
    
    % check if exponent matrix is the identity matrix
    if sum(sum(abs(E-diag(diag(E))))) == 0
        res = true;
    end

end

function res = aux_isPolytope(pZ,tol)

    res = true;
    
    % check if variable appears with exponent greater than 1
    if any(any(pZ.E >= 2))
       res = false;
       return;
    end

    % compute vertices
    [V,I] = aux_polyVertices(pZ);

    % loop over all vertices
    for i = 1:size(V,2)
       
        % compute vectors of normal cone
        N = aux_normalCone(pZ,I{i},tol);
        
        % loop over all vertices
        for j = 1:size(V,2)
           if i ~= j
               if ~aux_inCone(N,V(:,i),V(:,j),tol)
                   res = false;
                   return;
               end
           end
        end      
    end

end

function res = aux_inCone(N,v,v_,tol)
% checks if a vertex is inside the normal cone

    % construct constraints and objective function
    [n,m] = size(N);

    A1 = [N,-eye(n)];
    b1 = v_-v;
    A2 = [-N,-eye(n)];
    b2 = v-v_;
    
    lb = zeros(n+m,1);
    
    f = [zeros(1,m),ones(1,n)];
    
    % solve linear program
    persistent options
    if isempty(options)
        options = optimoptions('linprog','Display','off');
    end
    
    x = linprog(f,[A1;A2],[b1;b2],[],[],lb,[],options);
    
    % check if problem is solvable
    res = all(x(m+1:end) < tol);

end

function Nall = aux_normalCone(pZ,Ilist,tol)
% compute the normal cone for a vertex
    
    % loop over all factor combinations resulting in the same vertex
    Nall = [];
    
    for l = 1:size(Ilist,2)
        
        I = Ilist(:,l);
        N = zeros(size(pZ.G,1),length(I));

        % loop over all factors
        for i = 1:length(I)

            ind = find(pZ.E(i,:) > 0);

            % loop over all generators
            for j = 1:length(ind)
               k = ind(j);
               N(:,i) = N(:,i) + pZ.G(:,k) * prod(I.^pZ.E(:,k))/I(i)^pZ.E(i,k);
            end
        end
    
        N = N.*(-sign(I)');
        N = N(:,sum(abs(N),1) > tol);
        Nall = [Nall,N];
    end
    
end

function [V,Ilist] = aux_polyVertices(pZ)
% compute the polytope vertices and store the corresponding factor values

    % determine all potential vertices
    p = size(pZ.E,1);
    n = size(pZ.G,1);
    
    I = vertices(interval(-ones(p,1),ones(p,1)));
    
    V = zeros(n,size(I,2));
    
    for i = 1:size(I,2)
       V(:,i) = pZ.c;
       for j = 1:size(pZ.G,2)
           V(:,i) = V(:,i) + pZ.G(:,j) * prod(I(:,i).^pZ.E(:,j));
       end
    end
    
    % remove redundant points
    [V,ind] = sortrows(V');
    
    V = V';
    I = I(:,ind);

    V_ = zeros(size(V));
    Ilist = cell(size(V,2),1);
    
    V_(:,1) = V(:,1);
    Itemp = I(:,1);
    counter = 1;

    for i = 2:size(V,2)
        if ~all(abs(V(:,i)-V_(:,counter)) < 1e-10)
           counter = counter + 1;
           V_(:,counter) = V(:,i);
           Ilist{counter-1} = Itemp;
           Itemp = I(:,i);
        else
           Itemp = [Itemp,I(:,i)];
        end
    end
    
    Ilist{counter} = Itemp;

    V = V_(:,1:counter);
    Ilist = Ilist(1:counter);

    % determine vertices with the n-dimensional convex hull
    ind = convhulln(V');
    ind = unique(squeeze(ind));
    
    V = V(:,ind);
    Ilist = Ilist(ind);

end

% ------------------------------ END OF CODE ------------------------------
