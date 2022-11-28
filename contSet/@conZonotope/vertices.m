function V = vertices(cZ)
% vertices - calculate the vertices of a constrained zonotope object
%
% Syntax:  
%    V = vertices(cZ)
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
%    plotZono(cZ);
%    plot(V(1,:),V(2,:),'.k','MarkerSize',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_vertices('conZonotope',cZ);

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    V = vars{1}; return
end


% First remove redundant constraints that will not overapproximate the
% original constrained zonotope. Then we can check whether or not the
% resulting set is, in fact, a zonotope
cZ = reduceConstraints(cZ);

if ~isempty(cZ.A)
    
    % Calculate potential vertices of the constrained zonotope (vertices + 
    % points inside the set)
    V = potVertices(cZ);

    % Compute the convex hull to eliminate points located in the interior of
    % the constrained zonotope
    n = size(V,1);

    if size(V,2) > n+1    % set is full-dimensional
        ind = convhulln(V');
        ind = unique(ind,'stable');
        V = V(:,ind);
    end
    
else
    
   % no constraints -> call zonotope/vertices
   V = vertices(zonotope(cZ.Z));

end

%------------- END OF CODE --------------