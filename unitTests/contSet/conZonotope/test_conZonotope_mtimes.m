function res = test_conZonotope_mtimes
% test_conZonotope_mtimes - unit test function for the multiplication of a
%                           numerical or interval matrix with a contrained
%                           zonotope object
%
% Syntax:  
%    res = test_conZonotope_mtimes
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
% Written:      03-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;


% TEST 1: Random Test (numerical matrix, zonotope) ------------------------

% express a zonotope as a constrained zontope and compare the results of
% the mtimes operation for linear zonotopes and constrained zonotopes

for k = 1:5
   
    % generate a random zonotope object
    G = rand(2,10) - 0.5 * ones(2,10);
    zono = zonotope(G);
    cZono = conZonotope(G,[],[]);
    
    % generate a random numerical matrix
    mat = rand(2);
    
    % multiplication of matrix and zonotope
    zonoRes = mat * zono;
    cZonoRes = mat * cZono;
    
    % check if the results are equal
    
    if max(max(abs(zonoRes.Z - cZonoRes.Z))) > 1e-18
       error('Test 1 failed!');
    end   
end




% TEST 1: Random Test (interval matrix, zonotope) -------------------------

% express a zonotope as a constrained zontope and compare the results of
% the mtimes operation for linear zonotopes and constrained zonotopes

for k = 1:5
   
    % generate a random zonotope object
    G = rand(2,10) - 0.5 * ones(2,10);
    zono = zonotope(G);
    cZono = conZonotope(G,[],[]);
    
    % generate a random interval matrix
    m = rand(2) - 0.5*ones(2);
    r = rand(2);
    I = interval(m-r,m+r);
    
    % multiplication of interval matrix and zonotope
    zonoRes = I * zono;
    cZonoRes = I * cZono;
    
    % check if the results are equal
    
    if max(max(abs(zonoRes.Z - cZonoRes.Z))) > 1e-18
       error('Test 2 failed!');
    end   
end




% TEST 3: Random Test (interval matrix 2D) --------------------------------

for k = 1:5
    
    % Generate random conZonotope object
    points = rand(2,10);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    poly = mptPolytope(V');
    cZono = conZonotope(poly);

    % generate random interval matrix
    m = rand(2) - 0.5*ones(2);
    r = rand(2);
    I = interval(m-r,m+r);

    % Multiplication of interval matrix with constrained zonotope
    cZonoRes = I * cZono;

    % calculate extrem matrices + random matrices from the interval matrix
    ind = combinator(2,4);
    temp = cell(2,1);
    temp{1} = infimum(I);
    temp{2} = supremum(I);
    N = 100;

    matrices = cell(size(ind,1)+N);

    for i = 1:size(ind,1)
        matrices{i} = [temp{ind(i,1)}(1,1) temp{ind(i,2)}(1,2); ...
                       temp{ind(i,3)}(2,1) temp{ind(i,4)}(2,2)];
    end

    for i = 1:N
        temp = (rand(2) - 0.5*ones(2))*2;
        matrices{size(ind,1)+i} = m + temp.*r;
    end


    % calculate points that have to be located inside the resuling conZonotope
    N = size(V,2) * length(matrices);
    points = zeros(2,N);

    counter = 1;

    for i = 1:size(V,2)
        for j = 1:length(matrices)
            points(:,counter) = matrices{j} * V(:,i);
            counter = counter + 1;
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
%     plot(cZonoRes,[1,2],'r');

    % check if all points are located inside the resulting conZonotope
    Tol = 1e-14;

    for i = 1:size(points,2)
       if any(A*points(:,i) - b > Tol)
           error('Test 3 failed!');
       end
    end
end


res = 1;