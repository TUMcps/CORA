function res = test_conZonotope_enclose
% test_conZonotope_enclose - unit test function for the caclulation of the
%                            convex hull of two constrained zonotope 
%                            objects
%
% Syntax:  
%    res = test_conZonotope_enclose
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      28-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;


% TEST 1: Random Test (linear transformation 2D) --------------------------

for h = 1:5
    
    % generate random conZonotope object
    points = rand(2,10);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    poly = mptPolytope(V');
    cZono = conZonotope(poly);

    % generate random transformation matrix
    T = rand(2) - 0.5*ones(2);
    t = 5*(rand(2,1) - 0.5*ones(2,1));

    % transform the constrained zonotope
    cZono2 = T*cZono + t;

    % calculate convex hull
    cZonoRes = enclose(cZono,cZono2);

    % calculate points that have to be located inside the resuling conZonotope
    temp = vertices(cZono2);
    V2 = get(temp,'V');

    N = size(V,2) * size(V2,2) * 10;
    points = zeros(2,N);

    counter = 1;

    for i = 1:size(V,2)
        for j = 1:size(V2,2)
            for k = 0:0.1:1
                points(:,counter) = k * V(:,i) + (1-k) * V2(:,j);
                counter = counter + 1;
            end
        end
    end

    % convert the resulting conZonotope to a mptPolytope (to easily check if
    % a point is located inside the conZonotope)
    poly = mptPolytope(cZonoRes);

    % extract inequality constraints
    temp = get(poly,'P');
    A = temp.A;
    b = temp.b;

%     % visualize result
%     plot(points(1,:),points(2,:),'.k');
%     hold on
%     plot(cZono,[1,2],'b');
%     plot(cZono2,[1,2],'g');
%     plot(cZonoRes,[1,2],'r');

    % check if all points are located inside the resulting conZonotope
    Tol = 1e-14;

    for i = 1:size(points,2)
       if any(A*points(:,i) - b > Tol)
           error('Test 1 failed!');
       end
    end
end





% TEST 2: Random Test (zonotope addition 2D) ------------------------------

for h = 1:5
    
    % generate random conZonotope object
    points = rand(2,10);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    poly = mptPolytope(V');
    cZono = conZonotope(poly);

    % generate random transformation matrix
    T = rand(2) - 0.5*ones(2);
    t = 5*(rand(2,1) - 0.5*ones(2,1));
    
    % generate a random zonotope object
    zono = zonotope(rand(2,4)-0.5*ones(2,4));

    % transform the constrained zonotope and add the zonotope
    cZono2 = T*cZono + t + zono;

    % calculate convex hull
    cZonoRes = enclose(cZono,cZono2);

    % calculate points that have to be located inside the resuling conZonotope
    temp = vertices(cZono2);
    V2 = get(temp,'V');

    N = size(V,2) * size(V2,2) * 10;
    points = zeros(2,N);

    counter = 1;

    for i = 1:size(V,2)
        for j = 1:size(V2,2)
            for k = 0:0.1:1
                points(:,counter) = k * V(:,i) + (1-k) * V2(:,j);
                counter = counter + 1;
            end
        end
    end

    % convert the resulting conZonotope to a mptPolytope (to easily check if
    % a point is located inside the conZonotope)
    poly = mptPolytope(cZonoRes);

    % extract inequality constraints
    temp = get(poly,'P');
    A = temp.A;
    b = temp.b;

%     % visualize result
%     plot(points(1,:),points(2,:),'.k');
%     hold on
%     plot(cZono,[1,2],'b');
%     plot(cZono2,[1,2],'g');
%     plot(cZonoRes,[1,2],'r');

    % check if all points are located inside the resulting conZonotope
    Tol = 1e-14;

    for i = 1:size(points,2)
       if any(A*points(:,i) - b > Tol)
           error('Test 1 failed!');
       end
    end
end

res = 1;

