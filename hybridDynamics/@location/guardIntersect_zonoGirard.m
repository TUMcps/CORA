function R = guardIntersect_zonoGirard(loc,R,guard,options)
% guardIntersect_zonoGirard - implementation of the zonotope-hyperplane
%    intersection approach described in [1]
%
% Syntax:
%    R = guardIntersect_zonoGirard(loc,R,guard,options)
%
% Inputs:
%    loc - location object
%    R - list of intersections between the reachable set and the guard
%    guard - guard set (class: constrained hyperplane)
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - set enclosing the guard intersection
%
% References: 
%   [1] A. Girard et al. "Zonotope/Hyperplane Intersection for Hybrid 
%       Systems Reachablity Analysis"
%   [2] M. Althoff et al. "Zonotope bundles for the efficient computation 
%       of reachable sets", 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Stefan Liu, Niklas Kochdumper
% Written:       19-December-2016 
% Last update:   18-May-2018 (NK, integration into CORA)
%                19-December-2019 (NK, restructured the code)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check if guard set is a constrained hyperplane
    if ~isa(guard,'conHyperplane')
        throw(CORAerror('CORA:specialError',...
            "The method 'zonoGirard' only supports guards given as conHyperplane objects!")); 
    end

    % calc. orthogonal basis with the methods described in Sec. V.A in [2]
    B = calcBasis(loc,R,guard,options);
    
    % construct polytope from guard set inequality constraints C*x <= d
    if ~isempty(guard.C)
        poly = polytope(guard.C,guard.d);
    else
        poly = [];
    end
    
    % loop over all calculated basis 
    Z = cell(length(B),1);
    
    for i = 1:length(B)
       
        % loop over all reachable sets
        for j = 1:length(R)
           
            % interval enclosure in the transformed space according to [1]
            intTemp = aux_enclosingInterval(guard,B{i},R{j});
            
            % unite all intervals
            if j == 1
                I = intTemp;
            else
                I = I | intTemp;
            end
        end
        
        % set for one basis is empty -> overall set is emtpy
        if representsa_(I,'emptySet',eps)
            R = [];
            return;
        end
        
        % transform back to original space
        Z{i} = B{i}*zonotope(I);
        
        % remove parts outside the guard sets inequality constraints
        if ~representsa_(poly,"emptySet",eps)
            Z{i} = aux_tightenSet(Z{i},poly);
        end
    end
    
    % construct the enclosing zonotope bundle object
    Z = Z(~cellfun('isempty',Z));
    
    if isempty(Z)
        R = []; 
    elseif length(Z) == 1
        R = Z{1}; 
    else
        R = zonoBundle(Z); 
    end 
end


% Auxiliary functions -----------------------------------------------------

function I = aux_enclosingInterval(guard,B,Z)
% Implementation of Algorithm 2 in reference paper [1]

    % enclose the set with a zonotope
    if ~isa(Z,'zonotope') && ~isa(Z,'zonoBundle')
        Z = zonotope(Z);
    end

    % get hyperplane normal vector and offset
    n = guard.a';
    gamma = guard.b;

    % initialization
    lb = ones(length(n),1)*-inf;
    ub = ones(length(n),1)*inf;

    % loop over all basis vectors l 
    for i=1:length(B(1,:))

        if isa(Z,'zonoBundle')

            % loop over all parallel sets
            for k = 1:length(Z.Z)

                % Generate two-dimensional zonogon
                ZZ = [Z.Z{k}.c,Z.Z{k}.G];
                SZ = [(ZZ'*n)';(ZZ'*B(:,i))'];
                SZ(abs(SZ) < eps) = 0;
                Z_2D = zonotope(SZ);

                % Interval of intersection
                [m,M] = aux_bound_intersect_2D(Z_2D,gamma);

                lb(i) = max(m,lb(i));
                ub(i) = min(M,ub(i));
            end

        else
            % Generate two-dimensional zonogon
            ZZ = [Z.c,Z.G];
            SZ = [(ZZ'*n)';(ZZ'*B(:,i))'];
            SZ(abs(SZ) < eps) = 0;
            Z_2D = zonotope(SZ);

            % Interval of intersection
            [m,M] = aux_bound_intersect_2D(Z_2D,gamma);

            lb(i) = m;
            ub(i) = M;
        end
    end

    % Convert to Zonotope
    [lb,ub] = aux_robustProjection(B,n,gamma,lb,ub);

    % Increase robustness by equaling out small differences
    test = ub >= lb;

    if ~all(test)
        ind = find(test == 0);
        ub(ind) = ub(ind) + eps * ones(length(ind),1);
    end

    if all(lb <= ub)
        I = interval(lb,ub);
    else
        I = interval.empty(length(lb)); 
    end
end


function [m,M] = aux_bound_intersect_2D(Z,L)
% Implementation of Algorithm 3 in [1]

    Z = compact_(Z,'zeros',eps);

    c = center(Z);
    g = generators(Z);

    r = length(g(1,:));
    gamma = L; % L = {x,y: x = gamma}

    
    % Lower Bound ---------------------------------------------------------
    
    P=c; % current Position

    for i=1:r
        if (g(2,i) < 0) || ((g(2,i) == 0) && g(1,i) < 0)
            g(:,i) = - g(:,i); % ensure all generators are pointing upward
        end
        P = P - g(:,i); % drives P toward the lowest vertex of Z
    end

    % Direction of 
    if P(1) < gamma
        dir = 1;
        G = aux_sort_trig(g,dir); % we should look right
    else
        dir = -1;
        G = aux_sort_trig(g,dir); % or left 
    end
    dir = -dir;
    % G_low = G; %backup G

    s = sum(2*G,2);

    m = aux_dichotomicSearch(P,G,s,gamma,dir,Z);

    
    % Upper bound ---------------------------------------------------------
    
    P=c; % current Position

    for i=1:r
        if (g(2,i) < 0) || ((g(2,i) == 0) && g(1,i) < 0)
            g(:,i) = - g(:,i); % ensure all generators are pointing upward
        end
        P = P + g(:,i); % drives P toward the uppest vertex of Z
    end

    if P(1) > gamma
        dir = 1;
        G = -1*aux_sort_trig(g,dir); % we should look left
    else
        dir = -1;
        G = -1*aux_sort_trig(g,dir); % or right 
    end

    s = sum(2*G,2);

    M = aux_dichotomicSearch(P,G,s,gamma,dir,Z);

    if (abs(m) < 1e-12) && (abs(M) < 1e-12) %prevent numerical error
        m=0;M=0;
    elseif abs(m-M) <1e-12
        M=m;
    end
end

function G = aux_sort_trig(g,dir)
% sort g according to angle of every vector

    theta = cart2pol(g(1,:),g(2,:));
    if dir == 1
        [~,I] = sort(theta);
    elseif dir == -1
        [~,I] = sort(theta,'descend');
    end
    G = g(:,I);
end

function [G1,G2] = aux_split_pivot(G,s,dir)
% split the set of generators at the pivot element (= angle) 

    cos_G = dir*G(1,:)./sqrt(G(1,:).^2 + G(2,:).^2);
    pivot = dir*s(1)/sqrt(s(1)^2 + s(2)^2);

    G1 = G(:,cos_G <= pivot);
    G2 = G(:,cos_G > pivot);

end

function [m,m_over,m_unde] = aux_dichotomicSearch(P_0,G_0,s_0,gamma,dir,Z,queue)
% search the intersection between a 2D-zonotope and a line

    if nargin<7
        queue = cell(0,3);
    else
        queue(1,:) = [];
    end

    P=P_0;
    G=G_0;
    s=s_0;

    counter = 0;
    max_counter = length(G(1,:));

    while length(G(1,:)) > 1 && counter < max_counter
        [G1,G2] = aux_split_pivot(G,s,dir);
        s1 = sum(2*G1,2);
        line_int_y = aux_lineIntersect2D(P,P+s1,gamma);
        
        if ~isempty(line_int_y) % exist intersection
            G = G1;
            s = s1;
        else
            % save data for search
            new_entry = cell(1,3);
            new_entry{1,1} = P;
            new_entry{1,2} = G1;
            new_entry{1,3} = s1;
            queue = [queue;new_entry];

            G = G2;
            s = s-s1;
            P = P + s1;
        end
        counter = counter+1;
    end % only one generator remains

    m = aux_lineIntersect2D(P,P+s,gamma);

    if isempty(m)
        [m,m_over,m_unde] = aux_dichotomicSearch(...
            queue{1,1},queue{1,2},queue{1,3},gamma,dir,Z,queue);
    else
        % overapproximation
        m_over = P(2);
        % underapproximation
        m_unde = P(2)+s(2);
    end
end

function y = aux_lineIntersect2D(p1,p2,gamma)
% calculate the intersection between a line segment and a vertical line

    % check if gamma is one of the points
    if abs(p1(1) - gamma) < eps 
       y = p1(2);
       return;
    elseif abs(p2(1) - gamma) < eps
       y = p2(2);
       return;
    end

    % check if line is vertical
    if p1(1) == p2(1)
       y = [];
       return;
    end
    
    % check if line is horizontal
    if p1(2) == p2(2)
       if gamma > min(p1(1),p2(1)) && gamma < max(p1(1),p2(1))
          y = p2(2); 
       else
          y = []; 
       end
       return;
    end

    % calculate parameter of line y = a*x + b
    a = (p1(2)-p2(2))/(p1(1)-p2(1));
    b = p1(2)-a*p1(1);
    
    % calculate intersection with vertical line x = gamma
    y = a*gamma + b;
    
    % check if the intersection is located inside the line segment
    if y > max(p1(2),p2(2)) || y < min(p1(2),p2(2))
       y = []; 
    end
end

function [lb,ub] = aux_robustProjection(D,n,gamma,lb,ub)
% If the basis D is orthogonal to the hyperplane, set the interval value in
% the corresponding dimension to the hyperplane parameter to increase the
% numerical stability

    % project to modified space
    n_ = D' * n;
    
    % check if the basis is orhogonal to the hyperplane (with tolerance)
    ind = find(abs(n_) >= eps);
    
    % hyperplane normal vector perpendicular to orthogonal basis 
    % -> set the corresponding interval dimension to the hyperplane value
    if length(ind) == 1
        val = gamma * n_(ind);
        lb(ind) = val;
        ub(ind) = val;
    end
end

function Z = aux_tightenSet(Z,P)
% remove the parts of the interval that are located outside the inequality
% constraints C*x <= d of the constrained hyperplane

    % check if the zonotope fullfills the constraints
    if ~isIntersecting_(P,Z,'approx')
        Z = [];
        return;
    end
    
    % convert the zonotope to a constrained zonotope
    cZ = conZonotope(Z);
    
    % intersect the constrained zonotope with the polytope
    cZ = and_(cZ,P,'exact');
    
    % tighten the domain of the zonotope factors
    try
        cZ = rescale(cZ,'iter');
    catch
        cZ = rescale(cZ,'optimize');
    end
    
    % extract the rescaled zonotope from the constrained zonotope
    Z = zonotope(cZ.c,cZ.G);

end

% ------------------------------ END OF CODE ------------------------------
