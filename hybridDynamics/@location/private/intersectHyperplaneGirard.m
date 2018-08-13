function res = intersectHyperplaneGirard( Z, H, D )
% intersectHyperplaneGirard - intersect a zonotope with a hyperplane 
%
% Syntax:  
%    res = intersectHyperplaneGirard( Z, H, D )
%
% Inputs:
%    Z - zonotope
%    H - Hyperplane
%    D - cell array containing different orthonormal basis
%
% Outputs:
%    res - cell array containing the intervals determined with the
%          different orthonormal basis
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---
%
% References: 
%   [1] A. Girard et al. "Zonotope/Hyperplane Intersection for Hybrid 
%       Systems Reachablity Analysis"

% Author: Stefan Liu, Niklas Kochdumper
% Written: 19-Dec-2016 
% Last update: 18-May-2018 (NK, integration into CORA)
% Last revision: ---

%------------- BEGIN CODE --------------

% Implementation of Algorithm 2 in reference paper [1]

    % extract object properties
    if isa(H,'constrainedHyperplane')
        n = get(H.h,'c');
        gamma = get(H.h,'d');
    else
        n = get(H,'c');
        gamma = get(H,'d');
    end

    % initialization
    res.intervals = cell(length(D),1);

    % loop over all different orthogonal basis
    for j = 1:length(D)

        % initialization
        lb = ones(length(n),1)*-inf;
        ub = ones(length(n),1)*inf;
        
        % loop over all basis vectors l 
        for i=1:length(D{j}(1,:))
            
            if isa(Z,'zonotopeBundle')
                
                % loop over all parallel sets
                for k = 1:length(Z.Z)
                    
                    % Generate two-dimensional zonogon
                    ZZ = get(Z.Z{k},'Z');
                    SZ = [(ZZ'*n)';(ZZ'*D{j}(:,i))'];
                    SZ(abs(SZ) < eps) = 0;
                    Z_2D = zonotope(SZ);

                    % Interval of intersection
                    [m,M] = bound_intersect_2D(Z_2D,gamma);

                    lb(i) = max(m,lb(i));
                    ub(i) = min(M,ub(i));
                end
                
            else
                % Generate two-dimensional zonogon
                ZZ = get(Z,'Z');
                SZ = [(ZZ'*n)';(ZZ'*D{j}(:,i))'];
                SZ(abs(SZ) < eps) = 0;
                Z_2D = zonotope(SZ);

                % Interval of intersection
                [m,M] = bound_intersect_2D(Z_2D,gamma);

                lb(i) = m;
                ub(i) = M;
            end
        end

        % Convert to Zonotope
        [lb,ub] = robustProjection(D{j},n,gamma,lb,ub);
        
        temp = interval(lb,ub);
        
        if isempty(temp) % interval for one basis empty -> intersection empty
            res = [];
            break;
        else
            res.intervals{j} = temp;
        end
        
        % Check if the zonotope fullfills the constraints
        if isa(H,'constrainedHyperplane') && ...
           ~intersectsFeasibleRegion(zonotope(res.intervals{j}),H.C*D{j},H.d)
            res = [];
            break;
        end
    end
end





% Auxiliary Functions -----------------------------------------------------

function [ m,M ] = bound_intersect_2D( Z,L )
% Implementation of Algorithm 3 in reference paper [1]

    Z = deleteZeros(Z);

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
        G = sort_trig(g,dir); % we should look right
    else
        dir = -1;
        G = sort_trig(g,dir); % or left 
    end
    dir = -dir;
    % G_low = G; %backup G

    s = sum(2*G,2);

    m = dichotomicSearch(P,G,s,gamma,dir,Z);

    
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
        G = -1*sort_trig(g,dir); % we should look left
    else
        dir = -1;
        G = -1*sort_trig(g,dir); % or right 
    end

    s = sum(2*G,2);

    M = dichotomicSearch(P,G,s,gamma,dir,Z);

    if (abs(m) < 1e-12) && (abs(M) < 1e-12) %prevent numerical error
        m=0;M=0;
    elseif abs(m-M) <1e-12
        M=m;
    end
    
    % Debug
%     plot(Z,[1,2],'r');
%     hold on
%     v = vertices(Z);
%     V = get(v,'V');
%     mini = min(V(2,:));
%     maxi = max(V(2,:));
%     plot([gamma,gamma],[mini,maxi],'b');
%     plot(gamma,m,'.k','MarkerSize',15);
%     plot(gamma,M,'.k','MarkerSize',15);

end

function G = sort_trig(g,dir)
% sort g according to angle of every vector

    theta = cart2pol(g(1,:),g(2,:));
    if dir == 1
        [~,I] = sort(theta);
    elseif dir == -1
        [~,I] = sort(theta,'descend');
    end
    G = g(:,I);
end

function [G1,G2] = split_pivot(G,s,dir)
% split the set of generators at the pivot element (= angle) 

    cos_G = dir*G(1,:)./sqrt(G(1,:).^2 + G(2,:).^2);
    pivot = dir*s(1)/sqrt(s(1)^2 + s(2)^2);

    G1 = G(:,cos_G <= pivot);
    G2 = G(:,cos_G > pivot);

end

function [m,m_over,m_unde] = dichotomicSearch(P_0,G_0,s_0,gamma,dir,Z,queue)
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
        [G1,G2] = split_pivot(G,s,dir);
        s1 = sum(2*G1,2);
%         [~,line_int_y] = polyxpoly([P(1),P(1)+s1(1)],[P(2),P(2)+s1(2)],[gamma, gamma],[P(2),P(2)+s1(2)]);
        line_int_y = lineIntersect2D(P,P+s1,gamma);
        
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

%     [~,m] = polyxpoly([P(1),P(1)+s(1)],[P(2),P(2)+s(2)],[gamma, gamma],[P(2),P(2)+s(2)]);
    m = lineIntersect2D(P,P+s,gamma);

    if isempty(m)
        [m,m_over,m_unde] = dichotomicSearch(queue{1,1},queue{1,2},queue{1,3},gamma,dir,Z,queue);
    else
        m_over = P(2);% overapproximation
        m_unde = P(2)+s(2);% underapproximation
    end
end

function y = lineIntersect2D(p1,p2,gamma)
% calculate the intersection between a line segment and a vertical line

    % check if gamma is one of the points
    if p1(1) == gamma 
       y = p1(2);
       return;
    elseif p2(1) == gamma
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

function res = intersectsFeasibleRegion(zono,C,d)
% checks if parts of the zonotope are located inside the feasible region
% defined by the constraints

    res = 1;

    for i = 1:size(C,1)        
        
       % project zonotope onto normal direction of the hyperplane from the 
       % constraint
       proj = C(i,:) * zono.Z(:,1) - sum(abs(C(i,:) * zono.Z(:,2:end)));
       
       % check if zonotope is fully located outside the feasible region
       if proj > d(i)
           res = 0;
           break;
       end
    end
end

function [lb,ub] = robustProjection(D,n,gamma,lb,ub)
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

%------------- END OF CODE --------------