function V = vertices(obj)
% vertices - calculate the vertices of a constrained zonotope object
%
% Syntax:  
%    V = vertices(obj)
%
% Inputs:
%    obj - c-zonotope object
%
% Outputs:
%    V - matrix storing the vertices (dimension: [n,p], with p vertices)
%
% Example: 
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%    cZono = conZonotope(Z,A,b);
%
%    % calculate vertices
%    V = vertices(cZono);
%
%    % plot the result
%    hold on
%    plotZono(cZono);
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

% First remove redundant constraints that will not overapproximate the
% original constrained zonotope. Then we can check whether or not the
% resulting set is, in fact, a zonotope
obj = reduceConstraints(obj);

if ~isempty(obj.A)
    
    % Calculate potential vertices of the constrained zonotope (vertices + 
    % points inside the set)
    V = potVertices(obj);

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
   V = vertices(zonotope(obj.Z));
end


%------------- END OF CODE --------------