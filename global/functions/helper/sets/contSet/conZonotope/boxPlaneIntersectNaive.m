function L = boxPlaneIntersectNaive(hyperbox,hyperplane,beta)
% boxPlaneIntersectNaive - naive method of determining the intersection
%     points by the algorithm called Naive (Alg.1) provided in [1] going
%     along every edge and checking for intersections
%
% Mathematics:
%    edge between two adjacent vertices v and u (both dim n)
%    given by: e_uv(t) = v + t*(u-v) with t from 0 to 1
%    when the edge intersects the hyperplane the components
%    of e_uv sum up to alpha, hence t can be determined:
%          alpha - sum(i=1,n) v_i
%    t* = ------------------------
%           sum(i=1,n) u_i - v_i
%    if t* is between 0 and 1, there is an intersection
%
% Extension:
%    vector beta is not included above, hence:
%    own extension of the formula to include beta:
%          alpha - sum(i=1,n) beta_i * v_i
%    t* = ---------------------------------
%          sum(i=1,n) beta_i * (u_i - v_i)
%    if t* is between 0 and 1, there is an intersection
%
% Pseudocode from [1]:
%    INPUT:  hyperbox B and the hyperplane defined by alpha
%    OUTPUT: L ... list of solutions for the intersection
%    L <- {}
%    compute all edges E(B)
%    for all e_uv elem E(B) do
%        compute t*
%        if 0 <= t* <= 1 then
%            L <- L u e_uv(t*)
%        end if
%    end for
%    return L
%
% Syntax:
%    L = boxPlaneIntersectNaive(hyperbox,hyperplane,beta)
%
% Inputs:
%    hyperbox   - (n,2) with lower(:,1)/upper(:,2) bounds
%    hyperplane - val defined by constraining value alpha
%    beta       - scaling vector for hyperplane
%
% Outputs:
%    L - vertices of the intersection
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] C. Lara, J. Flores, F. Calderon.
%       "On the Hyperbox - Hyperplane Intersection Problem"

% Authors:       Mark Wetzlinger
% Written:       08-October-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1 define list of vertices
    L = [];
%   check if inputs are valid
%   mainly the equation: Aeq * x == b
%   e.g. (1 -0.5)*(x y)' == 2
%   1st check: alpha has to be a scalar
    if size(hyperplane) == 1
        alpha = hyperplane;
    else
        throw(CORAerror('CORA:specialError',...
            'The hyperplane has to be defined by a scalar.'));
    end   
%   2nd check: beta has to be a line vector (not a matrix)    
    if size(beta,1) > 1
        throw(CORAerror('CORA:specialError','Beta has to be a vector.'));
    end
% 2 compute all edges of hyperbox B
    dim = size(hyperbox,1);
    edges = aux_getEdges(dim);
    edges = aux_scaleEdges(edges,hyperbox);
% 3 begin loop
    for i=1:size(edges,1)
% 4     compute t*
        sumN = 0;
        sumD = 0;
        for d=1:dim
            sumN = sumN + beta(d) * edges(i,d,1);
            sumD = sumD + beta(d) * (edges(i,d,2) - edges(i,d,1));
        end
        t_star = (hyperplane - sumN)/sumD;
% 5     check if t* between 0 and 1
        if t_star >= 0 && t_star <= 1
% 6         add point to the list
            intersectionPoint = edges(i,:,1) + t_star*(edges(i,:,2) - edges(i,:,1));
%           additional constraint: should not already be in the list
            equality = 0;
            for k=1:size(L,1)
                if isequal(intersectionPoint,L(k,:))
                    equality = 1;
                    break;
                end
            end
            if equality == 0
                L = cat(1, intersectionPoint, L);
            end
% 7     end if
        end
% 8 end loop
    end
% 9 return list
end


% Auxiliary functions -----------------------------------------------------

function edges = aux_getEdges(dim)
%% description:
%   function calculates all adjacent vertices
%   given the dimension of the hyperbox
%   edges are abstracted to 0 and 1 (lower and upper bound)
%   start of edge is at edges(:,:,1), end is at edges(:,:,2)

%% code:
    points = zeros(2^dim,dim);
    edges = zeros(dim*2^(dim-1),dim,2);
    % index is dim which will be iterated next
    index = 1;
    % count goes through every vertex
    count = 1;
    % for position index in points and edges
    posP = 1;
    posE = 0;
    while count < 2^dim
        for i=1:dim
            if points(count,i) ~= 1
                newpoint = points(count,:);
                newpoint(i) = 1;
                if i >= index
                    posP = posP + 1;
                    points(posP,:) = newpoint;
                end
                posE = posE + 1;
                edges(posE,:,1) = points(count,:);
                edges(posE,:,2) = newpoint;
            end
        end
        % find index of first 0 after 1
        count = count + 1;
        index = find(points(count,:)) + 1;
    end
end

function edges = aux_scaleEdges(edges, B)
%% description:
%   scale 'edges' parametrized with 0 and 1 (lower and upper bound)
%   to the extension of the given hyperbox 'B'

%% code:
    numEdges = size(edges,1);
%   loop for start and end of every edge
    for i=1:2
%       loop for every point
        for j=1:numEdges
%           loop for every dimension
            for d=1:size(B,1)
                if edges(j,d,i) == 0
%                   0 -> use lower bound
                    edges(j,d,i) = B(d,1);
                else
%                   1 -> use upper bound
                    edges(j,d,i) = B(d,2);
                end
            end
        end 
    end
end

% ------------------------------ END OF CODE ------------------------------
