function V = vertices_(cZ,method,numDir,varargin)
% vertices_ - calculate the vertices of a constrained zonotope object
%
% Syntax:
%    V = vertices_(cZ)
%
% Inputs:
%    cZ - conZonotope object
%    method - 'default', 'template'
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
%                15-October-2024 (TL, integrated 'template' from plot function)
% Last revision: 27-March-2023 (MW, rename vertices_)

% ------------------------------ BEGIN CODE -------------------------------

% First, simplify the constrained zonotope as much as possible
cZ = compact(cZ, 'zeros', 1e-10);
% Check if generator matrix is empty
if isempty(cZ.G)
    V = cZ.c;
    return
end

% check if 2-dimensional
if dim(cZ) == 2
    % Vertices of a 2D cZ can be computed efficiently using support functions

    % compute vertices in projected dimensions
    V = projVertices(cZ,[1,2]);
    return
end

% First remove redundant constraints that will not overapproximate the
% original constrained zonotope. Then we can check whether or not the
% resulting set is, in fact, a zonotope
cZ = reduceConstraints(cZ);

% check if any constraint is left
if isempty(cZ.A)
   % no constraints -> call zonotope/vertices
   V = vertices(zonotope(cZ.c,cZ.G));
   return
end

% choose method
if strcmp(method,'template')
    V = aux_verticesTemplate(cZ,numDir);
    return
end

% default vertex computation
V = aux_verticesDefault(cZ);

end


% Auxiliary functions -----------------------------------------------------

function V = aux_verticesDefault(cZ)
    % Calculate potential vertices of the constrained zonotope
    % (vertices + points inside the set)
    V = priv_potVertices(cZ);
    
    % Compute the convex hull to eliminate points located in the interior of
    % the constrained zonotope
    n = size(V,1);
    
    % set is full-dimensional
    if rank(V,1e-6) > n+1
        ind = convhulln(V');
        ind = unique(ind,'stable');
        V = V(:,ind);
    end
end

function V = aux_verticesTemplate(cZ,numDir)
    % only implemented for 2d and 3d
    if ~any(dim(cZ) == [2,3])
        throw(CORAerror("CORA:notSupported",'Vertex computation using ''template'' is only supported for 2- and 3-dimensional constrained zonotopes.'))
    end

    % select directions for template polyhedron
    if dim(cZ) == 2
        % 2d
        angles = linspace(0,360,numDir+1);
        angles = angles(1:end-1);
        angles = deg2rad(angles);

        C = zeros(2,length(angles));

        for i = 1:length(angles)
            C(:,i) = [cos(angles(i));sin(angles(i))];
        end
    else
        % 3d
        N = ceil(sqrt(numDir));
        theta = 2 * pi * linspace(0,1,N);
        phi = acos(1 - 2 * linspace(0,1,N));
        [phi,theta] = meshgrid(phi,theta);
        C = zeros(3,numel(phi));
        cnt = 1;

        for i = 1:size(phi,1)
            for j = 1:size(phi,2)
                C(1,cnt) = sin(phi(i,j)) .* cos(theta(i,j));
                C(2,cnt) = sin(phi(i,j)) .* sin(theta(i,j));
                C(3,cnt) = cos(phi(i,j));
                cnt = cnt + 1;
            end
        end
    end

    % calculate the upper bounds along the directions
    d = zeros(size(C,2),1);
    for i = 1:length(d)
        d(i) = supportFunc_(cZ,C(:,i),'upper');
    end

    % compute intersection of neighboring constraints in 2D case
    if dim(cZ) == 2

        % loop over all pairs constraints to compute vertices
        V = zeros(2,numDir);
        for i=1:numDir
            if i==numDir
                % last constraint with first constraint
                V(:,i) = C(:,[i,1])' \ d([i,1]);
            else
                V(:,i) = C(:,[i,i+1])' \ d(i:i+1);
            end
        end
        % remove duplicates using relative tolerance
        [V,IA] = uniquetol(V',1e-3,'ByRows',true);
        % re-order
        [~,order] = mink(IA,length(IA));
        V = V(order,:)';

    else

        % construct template polyhedron
        poly = polytope(C',d);

        % ######

    end

end

% ------------------------------ END OF CODE ------------------------------
