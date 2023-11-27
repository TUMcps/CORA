function V = removeCollinearVertices2D(V,varargin)
% removeCollinearVertices2D - removes collinear points up to a user-defined
%    tolerance by checking the rank of two vectors computed from two
%    subsequent pairs of vertices in the list
%
% Syntax:
%    V = removeCollinearVertices2D(V)
%    V = removeCollinearVertices2D(V,tol)
%
% Inputs:
%    V - vertices
%    tol - (optional) tolerance
%
% Outputs:
%    V - vertices
%
% Example:
%    V = 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default tolerance
tol = setDefaultValues({0}, varargin);

% number of vertices
numVert = size(V,2);

% too few vertices
if size(V,2) <= 2
    return
end

% order vertices by angle
meanV = mean(V,2);
angle = atan2(V(2,:)-meanV(2),V(1,:)-meanV(1));
[~,angleIdx] = sort(angle);
V = V(:,angleIdx);

% compute vectors (note: wrap around at the end)
vectors = [V, V(:,1)] - [V(:,end), V];

% indices for which vertices are kept
idxKept = true(1,numVert);
cnt = 0;

% check rank
for i = 1:numVert
    if rank(vectors(:,i:i+1), tol) < 2
        % remove i-th vertex
        idxKept(i) = false;

        % increment counter of removed vertices
        cnt = cnt + 1;

        % if only two vertices remain, abort
        if (numVert - cnt) == 2
            break
        end
    end
end

% remove collinear vertices (if any)
if ~all(idxKept)
    V = V(:,idxKept);
end

% ------------------------------ END OF CODE ------------------------------
