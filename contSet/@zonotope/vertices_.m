function V = vertices_(Z,alg,varargin)
% vertices_ - returns potential vertices of a zonotope
%    WARNING: Do not use this function for high-order zonotopes as the
%    computational complexity grows exponentially!
%
% Syntax:
%    V = vertices_(Z)
%    V = vertices_(Z,alg)
%
% Inputs:
%    Z - zonotope object
%    alg - algorithm used
%           - 'convHull' (default)
%           - 'iterate'
%           - 'polytope'
%
% Outputs:
%    V - matrix
%
% Example:
%    Z = zonotope([1;-1],[1 3 -2 1 0; 0 2 1 -2 1]);
%    V = vertices(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices, interval, polytope

% Authors:       Matthias Althoff, Niklas Kochdumper, Daniel Heß
% Written:       14-September-2006 
% Last update:   30-January-2008
%                23-June-2009
%                24-June-2010
%                11-July-2012
%                28-June-2016 (DH, implemented fast 2D computation)
%                30-July-2016
%                28-October-2019 (NK, added new algorithm)
%                02-June-2023 (TL, convHull with degenerate zonotopes)
%                11-October-2024 (TL, moved fast 2D method here)
% Last revision: 27-March-2023 (MW, rename vertices_)

% ------------------------------ BEGIN CODE -------------------------------

% different cases for different dimensions
n = dim(Z);

if n == 1
    % compute the two vertices for one-dimensional case
    c = Z.c;
    temp = sum(abs(Z.G),2);
    V = [c - temp,c + temp];

elseif n == 2
    % use fast method
    V = aux_vertices2Dfast(Z);

else
    % apply the selected algorithm
    if strcmp(alg,'iterate')
        V = aux_verticesIterate(Z);
    elseif strcmp(alg,'polytope')
        V = aux_verticesPolytope(Z);
    else
        V = aux_verticesConvHull(Z);
    end
end

% remove duplicates
V = unique(V','rows','stable')';
    
end


% Auxiliary functions -----------------------------------------------------

function V = aux_vertices2Dfast(Z)
    % use fast method for 2D zonotopes
    % (this function was previously in zonotope/polygon)
    
    % empty case
    if representsa_(Z,'emptySet',eps)
        V = zeros(dim(Z),0); return
    end
    
    % delete zero generators
    Z = compact_(Z,'zeros',eps);
    
    % obtain center and generator matrix
    c = Z.c;
    G = Z.G;
    
    % obtain number of generators
    nrGens = size(G,2);
    
    % obtain size of enclosing intervalhull of first two dimensions
    xmax = sum(abs(G(1,:)));
    ymax = sum(abs(G(2,:)));
     
    % Z with normalized direction: All generators pointing "up"
    Gnorm = G;
    Gnorm(:,G(2,:)<0)=G(:,G(2,:)<0)*-1;
    
    %compute angles
    angles = atan2(Gnorm(2,:),Gnorm(1,:));
    angles(angles<0) = angles(angles<0) +2*pi;%handle numerical imprecision/deficiency in atan2, wraparound is not at pi?!?
    
    % assert(~any(angles>pi));
    
    %sort all generators by their angle
    [~,IX] = sort(angles,'ascend');
    
    %cumsum the generators in order of angle
    V = zeros(2,nrGens+1);
    for i = 1:nrGens
        V(:,i+1) = V(:,i) + 2*Gnorm(:,IX(i));
    end
    
    V(1,:) = V(1,:) + xmax - max(V(1,:));
    V(2,:) = V(2,:) - ymax;
    
    %flip/mirror upper half to get lower half of zonotope (point symmetry)            
    V = [V(1,:),V(1,end)+V(1,1)-V(1,2:end);...
        V(2,:),V(2,end)+V(2,1)-V(2,2:end)];
    
    %consider center
    V = c + V;

end

function V = aux_verticesPolytope(Z)

    P = polytope(Z);
    V = vertices(P);

end

function V = aux_verticesConvHull(Z)

    % first vertex is the center of the zonotope
    c = Z.c;
    G = Z.G;
    n = dim(Z);
    nrGens = size(G,2);

    % we project vertices to lower-dims if degenerate
    % starting with 1 points (center); always degenerate
    V = c;
    last_n = 0;

    % generate further potential vertices in the loop
    for i=1:nrGens

        % expand vertices with current generator
        g = G(:,i);
        V = [V-g,V+g];

        % remove inner points
        try
            if last_n < n
                % last iteration it was degenerate, 
                % assuming it still is; project to lower-dim

                % shift by center
                c = Z.c;
                
                % compute projection matrix via SVD
                [U,S] = svd(V-c,'vector');
                P = U(:,~withinTol(S,0));
                
                % project into low-dim space
                V_proj = P'*V;
    
                % check if still degenerate
                last_n = size(V_proj,1);
    
                if last_n == 1
                    % 1-dim can be computed quickly
                    [~,imin] = min(V_proj); 
                    [~,imax] = max(V_proj); 
                    indices = [imin,imax];
                else
                    % compute convex hull
                    K = convhulln(V_proj');
                    indices = unique(K);
                end
                
                % select corresponding points of original matrix
                V = V(:,indices);

            else
                % non-degenerate; normal conv hull computation
                K = convhulln(V');
                indices = unique(K);
                V = V(:,indices);
            end
            
        catch ME
            CORAwarning('CORA:contSet','Convex hull failed. Continuing...')
        end
    end
end

function V = aux_verticesIterate(Z)

    % delete aligned and all-zero generators
    Z = compact_(Z,'zeros',eps);
    Z = compact_(Z,'aligned',1e-3);

    % extract object data
    G = generators(Z);
    c = center(Z);
    n = size(G,1);
    
    % catch the case where the zonotope is not full-dimensional
    if size(G,2) < n
        V = aux_verticesIterateSVG(Z);
        return;
    end

    % compute vertices of the parallelotope
    vert = vertices_(interval(-ones(n,1),ones(n,1)));  
    V = c + G(:,1:n)*vert;
    
    % compute halfspaces of the parallelotope
    [poly,~,isDeg] = polytope(zonotope([c,G(:,1:n)]));
    if isDeg
        V = aux_verticesIterateSVG(Z);
        return;
    else
        A = poly.A;
    end
    
    % loop over all remaining generators 
    for i = n+1:size(G,2)
       
        % extract current generator
        g = G(:,i);
        
        % compute potential vertices
        V = [V+g,V-g];
        
        % compute new halfspaces
        if n == 2
            temp = ndimCross(g);
            temp = temp/norm(temp);
            Anew = [temp';-temp'];
        else
            comb = combinator(i-1,n-2,'c');
            Anew = zeros(2*size(comb,1),n);
            counter = 1;

            for j = 1:size(comb,1)
                temp = ndimCross([G(:,comb(j,:)),g]);
                temp = temp/norm(temp);
                Anew(counter,:) = temp';
                Anew(counter+1,:) = -temp';
                counter = counter + 2;
            end
        end
        
        A = [A;Anew];
        
        % compute halfspace offsets
        b = max(A*V,[],2);
        
        % remove redundant vertices
        temp = max(A*V-b,[],1); 
        nV = aux_numVertices(i,n);
        [~,ind] = sort(temp,'descend');
        V = V(:,ind(1:nV));
        
    end

end

function [res,suc] = aux_verticesIterateSVG(Z)
% compute vertices for the case that zonotope is not full-dimensional

    suc = true;
    res = [];

    % extract object data
    G = generators(Z);
    c = center(Z);
    n = dim(Z);

    % singular value decomposition
    [S,V,~] = svd(G);
    
    if size(V,2) < size(G,1)
       V = [V,zeros(n,n-size(V,2))]; 
    end

    % state space transformation
    Z_ = S'*[c,G];

    % remove dimensions with all zeros
    ind = find(diag(V) <= 1e-12);
    ind_ = setdiff(1:size(V,1),ind);

    if ~isempty(ind)
        % compute vertices in transformed space
        V = vertices(zonotope(Z_(ind_,:)));

        % transform back to original space
        res = S*[V;zeros(length(ind),size(V,2))];
    else
        suc = false;
    end 
end

function res = aux_numVertices(m,n)
% compute number of zonotope vertices

    res = 0;
    for i = 0:n-1
        res = res + nchoosek(m-1,i); 
    end
    res = 2*res;
end

% ------------------------------ END OF CODE ------------------------------
