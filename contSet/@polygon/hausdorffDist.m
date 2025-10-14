function val = hausdorffDist(pgon,S)
% hausdorffDist - Calculates the Hausdorff distance between two polygons or 
%    a polygon and a point
%
% Syntax:
%    val = hausdorffDist(pgon,S)
%
% Inputs:
%    pgon - polygon object
%    S - contSet object of single point
%
% Outputs:
%    val - Hausdorff distance
%
% Examples:
%    pgon1 = polygon.generateRandom()
%    pgon2 = polygon.generateRandom()
%   
%    val = hausdorffDist(pgon1,pgon2)
%
%    figure; hold on;
%    plot(pgon1);
%    plot(pgon2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/hausdorffDist

% Authors:       Niklas Kochdumper
% Written:       28-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init result
val = 0;

% read out polytope object
[pgon,S] = findClassArg(pgon,S,'polygon');

% check dimensions
equalDimCheck(pgon,S);

% different cases for different types of sets
if isnumeric(S)

    % catch special cases
    if representsa_(pgon,'emptySet')
        if isempty(S)
            val = 0; return
        else
            val = Inf; return
        end
    end

    % compute Hausdorff distance for each point    
    for i = 1:size(S,2)
        if contains(pgon,S(:,i))
            val_ = 0;
        else
            val_ = aux_distPolygonPoint(pgon,S(:,i));
        end
        val = max(val,val_);
    end
    
elseif isa(S,'polygon') || isa(S,'interval') || ...
       isa(S,'zonotope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle') || isa(S,'polytope')

    % convert set to polygon
    S = polygon(S);

    % loop over all vertices of the first polygon
    V = vertices(pgon);

    for i = 1:size(V,2)
        val_ = aux_distPolygonPoint(S,V(:,i));
        val = max(val,val_);
    end

    % loop over all vertices of the second polygon
    V = vertices(S);

    for i = 1:size(V,2)
        val_ = aux_distPolygonPoint(pgon,V(:,i));
        val = max(val,val_);
    end
   
else
    % throw error for given arguments
    throw(CORAerror('CORA:noops',pgon,S));
end

end


% Auxiliary functions -----------------------------------------------------

function val = aux_distPolygonPoint(P,x)
% compute Hausdorff distance between a polygon and a single point

    % loop over the polygon boundary
    V = vertices(P); V = [V,V(:,1)]; val = inf;

    for i = 1:size(V,2)-1
        val_ = aux_distPointLineSegment(V(:,i),V(:,i+1),x);
        val = min(val,val_);
    end
end

function val = aux_distPointLineSegment(l1,l2,x)
% comptute the minimum distance between a line segment formed by the points
% l1 and l2 and the point x

    % minimum distance to points l1 and l2
    val = min(norm(l1-x),norm(l2-x));

    % compute hyperplane representation of the line segment
    d = l2 - l1; n = [d(2);-d(1)]; c = n'*l1;

    % project point x onto the hyperplane
    p = x + n*(c - n'*x)/(n'*n);

    % check if projected point is on the line segment
    if d'*p >= d'*l1 && d'*p <= d'*l2
        val = min(val,norm(p-x));
    end
end

% ------------------------------ END OF CODE ------------------------------
