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
%    V - vertices object
%
% Example: 
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%    cZono = conZonotope(Z,A,b);
%
%    % calculate vertices
%    v = vertices(cZono);
%    V = get(v,'V');
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
    
if ~isempty(obj.A)
    
    % Calculate potential vertices of the constrained zonotope (vertices + 
    % points inside the set)
    z = potVertices(obj);

    % Compute the convex hull to eliminate points located in the interior of
    % the constrained zonotope
    dim = size(z,1);

    if size(z,2) > dim+1    % set is full-dimensional
        ind = convhulln(z');
        ind = unique(ind,'stable');
        z = z(:,ind);
    end

    % Construct vertices object
    V = vertices(z);
    
else
    
   % no constraints -> call zonotope/vertices
   V = vertices@zonotope(obj);
end


%------------- END OF CODE --------------