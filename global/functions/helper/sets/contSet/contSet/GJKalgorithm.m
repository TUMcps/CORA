function res = GJKalgorithm(S1,S2,varargin)
% GJKalgorithm - check for intersection between two convex sets using the
%    Gilbert-Johnson-Keerthi algorithm [1]
%
% Syntax:
%    res = GJKalgorithm(S1,S2)
%    res = GJKalgorithm(S1,S2,tol)
%
% Inputs:
%    S1 - contSet object representing a convex set
%    S2 - contSet object representing a convex set
%    tol - tolerance
%
% Outputs:
%    res - true is sets intersect, false if not
%
% Example:
%    S1 = polytope([0.5 2 4; 1 3 1]);
%    S2 = polytope([2.5 3 6 6; 2 3.5 3.5 1]);
% 
%    GJKalgorithm(S1,S2)
%
%    figure; hold on; box on;
%    plot(S1,[1,2],'r');
%    plot(S2,[1,2],'b');
%
% References:
%    [1] E. Gilbert, D. Johnson, and S. Keerthi. "A fast procedure for 
%        computing the distance between complex objects in 
%        three-dimensional space", 1988
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isIntersecting

% Authors:       Niklas Kochdumper
% Written:       21-April-2022
% Last update:   24-April-2024 (MW, include tolerance, enforce convexity)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % ensure that dimensions match
    equalDimCheck(S1,S2);

    % set tolerance
    tol = setDefaultValues({1e-14},varargin);

    % ensure that both sets are convex
    if ~representsa_(S1,'convexSet',tol) || ~representsa_(S2,'convexSet',tol)
        throw(CORAerror('CORA:notSupported',['The GJK algorithm can only ' ...
            'be called for convex sets.']));
    end

    % read out dimension
    n = dim(S1);

    % initial direction: vector between centers or first basis vector
    try
        d = center(S2) - center(S1);
    catch
        d = zeros(n,1); d(1) = 1;
    end

    % compute initial support vector
    [~,s1] = supportFunc(S1,d);
    [~,s2] = supportFunc(S2,-d);
    v = s1 - s2;

    % check for immediate exit (vertex is origin)
    if all(withinTol(v,zeros(n,1),tol))
        res = true; return
    end
    
    % initialize simplex (in vertex representation) and next direction
    simplex = v;
    d = -v;
    
    % loop until intersection can be proven or disproven
    while true
       
        % compute support vector
        [~,s1] = supportFunc(S1,d);
        [~,s2] = supportFunc(S2,-d);
        v = s1 - s2;
        
        % check if origin crossed -> no intersection if not
        if d'*v < 0
            res = false; return 
        end
        
        % add current support vector to simplex
        simplex = [simplex, v];
        
        % determine simplex face that is closes to the origin
        [simplex,d,containsOrigin] = aux_nearestSimplex(simplex);
        
        if containsOrigin
            res = true; return;
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function [simplex,d,containsOrigin] = aux_nearestSimplex(simplex)
% determine new search direction that is orthogonal to the face to the 
% simplex (given as a polytope in vertex representation, in this case, a
% matrix) that is closest to the origin
% returns: simplex in vertex representation, next direction, exit flag

    containsOrigin = false;
    n = size(simplex,1);

    % check if the simplex is already full-dimensional
    if size(simplex,2) > size(simplex,1)
        
        % convert the simplex to halfspace representation
        P = polytope(simplex);
        
        % if the simplex contains the origin, the GJK algorithm has
        % determined that the sets intersect -> exit immediately
        if contains(P,zeros(n,1))
           containsOrigin = true; d = []; return; 
        end       
        
        % use different algorithms for computing the minimum distance 
        % between the simplex and the origin depending on the dimension
        if n > 3
            index = aux_minDistOriginQuadProg(P);
        else
            index = aux_minDistOriginRecursive(P);
        end
        
        % construct next simplex and next search direction
        d = P.A(index,:)';
        
        [~,ind] = sort(abs(P.A(index,:)*simplex - P.b(index)));
        simplex = simplex(:,ind(1:n));
        
    else
        
        % not full-dimensional -> determine new direction that points
        % towards the origin and is orthogonal to the current simplex
        s = [simplex, zeros(n,1)];
        orthBasis = gramSchmidt(s(:,2:end) - s(:,1:end-1));
        d = orthBasis(:,end);
    end
end

function index = aux_minDistOriginQuadProg(P)
% compute the minimum distance between a simplex (represented as a
% polytope) and the origin using quadratic programming

    % object properties
    A = P.A; b = P.b; ind = 1:size(A,1); n = dim(P);

    % check which halfspaces point toward the origin
    queue = find(b <= 0);
    
    if length(queue) == 1
       index = queue(1); return; 
    end

    % loop over all halspaces to find the one that is closest to origin
    dist = inf; options = optimoptions('quadprog','Display','none');

    for i = 1:length(queue)
        ind_ = setdiff(ind,queue(i));
        [~,dist_] = quadprog(eye(n),zeros(n,1),A(ind_,:),b(ind_), ...
                             A(queue(i),:),b(queue(i)),[],[],[],options);
        if dist_ < dist
           dist = dist_; index = queue(i); 
        end
    end
end

function [index,point] = aux_minDistOriginRecursive(P)
% compute the minimum distance between a simplex and the origin by
% recusively decomposing the simplex into lower dimensional simplices

    % initialization
    dist = inf;
    A = P.A; b = P.b; V = vertices(P); n = dim(P);
    
    % one-dimensional simplex
    if size(V,1) == 1
       [point,index] = min(V); return;
    end

    % loop over all polytope halfspaces
    nrCon = size(A,1);
    for i = 1:nrCon
        
        if b(i) <= 0
        
            % construct orthonormal basis
            B = gramSchmidt(A(i,:)');

            % construct n-1-dimensional simplex corresponding to halfspace
            [~,ind] = sort(abs(A(i,:)*V - b(i)));
            V_ = B'*V(:,ind(1:n));
            poly_ = polytope(V_(2:end,:));

            % point on the hyperplane that is closest to the origin
            point_ = b(i)/norm(A(i,:))^2 * A(i,:)';
            p = B' * point_;

            % check if closest point is in the simplex -> minimum distance 
            if ~contains(poly_,p(2:end))
                [~,p_] = aux_minDistOriginRecursive(poly_ - p(2:end)); 
                point_ = point_ + B * [0; p_];
            end

            % update minimum distance
            dist_ = norm(point_);

            if dist_ < dist
               dist = dist_; point = point_; index = i; 
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
