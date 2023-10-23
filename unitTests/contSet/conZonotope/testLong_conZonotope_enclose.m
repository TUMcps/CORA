function res = testLong_conZonotope_enclose
% testLong_conZonotope_enclose - unit test function for the
%   calculation of the convex hull of two constrained zonotope objects
%
% Syntax:
%    res = testLong_conZonotope_enclose
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
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

% Authors:       Niklas Kochdumper
% Written:       28-June-2018
% Last update:   28-April-2019 (MA, code shortened)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;


% TEST 1: Random Test (linear transformation 2D) --------------------------

for h = 1:5
    
    % generate random conZonotope object
    points = rand(2,10);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    P = polytope(V);
    cZono = conZonotope(P);

    % generate random transformation matrix
    T = rand(2) - 0.5*ones(2);
    t = 5*(rand(2,1) - 0.5*ones(2,1));

    % transform the constrained zonotope
    cZono2 = T*cZono + t;

    % calculate convex hull
    cZonoRes = enclose(cZono,cZono2);

    % calculate points that have to be located inside the resulting conZonotope
    V2 = vertices(cZono2);

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

    % convert the resulting conZonotope to a polytope (to easily check if
    % a point is located inside the conZonotope)
    P = polytope(cZonoRes);

    % extract inequality constraints
    A = P.A;
    b = P.b;

%     % visualize result
%     plot(points(1,:),points(2,:),'.k');
%     hold on
%     plot(cZono,[1,2],'b');
%     plot(cZono2,[1,2],'g');
%     plot(cZonoRes,[1,2],'r');

    % check if all points are located inside the resulting conZonotope
    for i = 1:size(points,2)
       if ~all(A*points(:,i) < b | withinTol(A*points(:,i),b))
          throw(CORAerror('CORA:testFailed'));
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

    P = polytope(V);
    cZono = conZonotope(P);

    % generate random transformation matrix
    T = rand(2) - 0.5*ones(2);
    t = 5*(rand(2,1) - 0.5*ones(2,1));
    
    % generate a random zonotope object
    zono = zonotope(rand(2,4)-0.5*ones(2,4));

    % transform the constrained zonotope and add the zonotope
    cZono2 = T*cZono + t + zono;

    % calculate convex hull
    cZonoRes = enclose(cZono,cZono2);

    % calculate points that have to be located inside the resulting conZonotope
    V2 = vertices(cZono2);

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

    % convert the resulting conZonotope to a polytope (to easily check if
    % a point is located inside the conZonotope)
    P = polytope(cZonoRes);

    % extract inequality constraints
    A = P.A;
    b = P.b;

%     % visualize result
%     plot(points(1,:),points(2,:),'.k');
%     hold on
%     plot(cZono,[1,2],'b');
%     plot(cZono2,[1,2],'g');
%     plot(cZonoRes,[1,2],'r');

    % check if all points are located inside the resulting conZonotope
    for i = 1:size(points,2)
       if ~all(A*points(:,i) < b | withinTol(A*points(:,i),b))
          throw(CORAerror('CORA:testFailed'));
       end
    end
end

% ------------------------------ END OF CODE ------------------------------
