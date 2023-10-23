function V = vertices_(cZ,varargin)
% vertices_ - calculate the vertices of a constrained zonotope object
%
% Syntax:
%    V = vertices_(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    V - matrix storing the vertices (dimension: [n,p], with p vertices)
%
% Example: 
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%    V = vertices(cZ);
% 
%    figure; hold on
%    plot(cZ);
%    plot(V(1,:),V(2,:),'.k','MarkerSize',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices

% Authors:       Niklas Kochdumper
% Written:       11-May-2018
% Last update:   25-April-2023 (TL, 2d support func computation)
% Last revision: 27-March-2023 (MW, rename vertices_)

% ------------------------------ BEGIN CODE -------------------------------

if dim(cZ) == 2
    % Vertices of a 2D cZ can be computed efficiently using support functions
        
    % compute vertices in projected dimensions
    V = projVertices(cZ,[1,2]);
else
        
    % First remove redundant constraints that will not overapproximate the
    % original constrained zonotope. Then we can check whether or not the
    % resulting set is, in fact, a zonotope
    cZ = reduceConstraints(cZ);
    
    if ~isempty(cZ.A)
        
        % Calculate potential vertices of the constrained zonotope
        % (vertices + points inside the set)
        V = potVertices(cZ);
    
        % Compute the convex hull to eliminate points located in the interior of
        % the constrained zonotope
        n = size(V,1);
    
        % set is full-dimensional
        if rank(V,1e-6) > n+1
            ind = convhulln(V');
            ind = unique(ind,'stable');
            V = V(:,ind);
        end
        
    else
        
       % no constraints -> call zonotope/vertices
       V = vertices(zonotope(cZ.c,cZ.G));
    
    end
end

% ------------------------------ END OF CODE ------------------------------
