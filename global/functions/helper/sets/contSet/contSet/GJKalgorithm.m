function res = GJKalgorithm(S1,S2)
% GJKalgorithm - check for intersection between two convex sets using the
%                Gilbert-Johnson-Keerthi algorithm
%
% Syntax:
%    res = GJKalgorithm(S1,S2)
%
% Inputs:
%    S1 - contSet object representing a convex set
%    S2 - contSet object representing a convex set
%
% Outputs:
%    res - true is sets intersect, false if not
%
% Example:
%    
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    n = dim(S1);

    % get initial direction
    try
        d = center(S2) - center(S1);
    catch
        d = zeros(n,1); d(1) = 1;
    end

    % compute initial support vector
    [~,s1] = supportFunc(S1,d);
    [~,s2] = supportFunc(S2,-d);
    v = s1 - s2;
    
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
% simplex that is closest to the origin

    containsOrigin = false;
    n = size(simplex,1);

    % check if the simplex is already full-dimensional
    if size(simplex,2) > size(simplex,1)
        
        % convert the simplex to halfspace representation
        poly = polytope(simplex);
        
        if contains(poly,zeros(n,1))
           containsOrigin = true; d = []; return; 
        end       
        
        % use different algorithms for computing the minimum distance 
        % between a simplex and the origin depending on the dimension
        if n > 3
            index = aux_minDistOriginQuadProg(poly);
        else
            index = aux_minDistOriginRecursive(poly);
        end
        
        % construct next simplex and next search direction
        d = poly.A(index,:)';
        
        [~,ind] = sort(abs(poly.A(index,:)*simplex - poly.b(index)));
        simplex = simplex(:,ind(1:n));
        
    else
        
        % not full-dimensional -> determine new direction that points
        % towards the origin and is orthogonal to the current simplex
        s = [simplex, zeros(n,1)];
        tmp = gramSchmidt(s(:,2:end) - s(:,1:end-1));
        d = tmp(:,end);
    end
end

function index = aux_minDistOriginQuadProg(poly)
% compute the minimum distance between a simplex and the origin using
% quadratic programming

    % object properties
    A = poly.A; b = poly.b; ind = 1:size(A,1); n = size(A,2);

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

function [index,point] = aux_minDistOriginRecursive(poly)
% compute the minimum distance between a simplex and the origin by
% recusively decomposing the simplex into lower dimensional simplices

    % initialization
    dist = inf;
    A = poly.A; b = poly.b; V = vertices(poly); n = size(A,2);
    
    % one-dimensional simplex
    if size(V,1) == 1
       [point,index] = min(V); return;
    end

    % loop over all polytope halfspaces
    for i = 1:size(A,1)
        
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
