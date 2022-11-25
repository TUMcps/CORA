function res = testLongDuration_polyZonotope_mptPolytope
% testLongDuration_polyZonotope_mptPolytope - unit test function for the
%    conversion of a mptPolytope to a polynomial zonotope
%
% Syntax:  
%    res = test_polyZonotope_mptPolytope
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

% Author:       Niklas Kochdumper
% Written:      30-October-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;

%% RANDOM TESTS

% TEST 2-dimensional (Polytope -> Polynomial Zonotope)

for j = 10:15
    
    % Generate random polytope vertices
    points = rand(2,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind(1:min(j,length(ind))));

    % Construct a mptPolytope object from the vertices
    poly = mptPolytope(V');

    % Convert to polynomial zonotope
    pZres = polyZonotope(poly);
    
    % check if the all polytope vertices are located inside the polynomial
    % zonotope
    suc = containsPointSet(pZres,V,[],30);
    
    if ~suc
       error('test_polyZonotope_mptPolytope: random test 2D (polytope -> polyZonotope) failed!'); 
    end   
end



% TEST 4-dimensional (Polytope -> Polynomial Zonotope)

for j = 10:12
    
    % Generate random polytope vertices
    points = rand(4,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind(1:min(j,length(ind))));

    % Construct a mptPolytope object from the vertices
    poly = mptPolytope(V');

    % Convert to polynomial zonotope
    pZres = polyZonotope(poly);
    
    % check if the all polytope vertices are located inside the polynomial
    % zonotope
    suc = containsPointSet(pZres,V,[],30);
    
    if ~suc
       error('test_polyZonotope_mptPolytope: random test 4D (polytope -> polyZonotope) failed!'); 
    end   
end



% TEST 2-dimensional (Polynomial Zonotope -> Polytope)

for j = 10:15

    % Generate random polytope vertices
    points = rand(2,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind(1:min(j,length(ind))));

    % Construct a mptPolytope object from the vertices
    poly = mptPolytope(V');

    % Convert to polynomial zonotope
    pZ = polyZonotope(poly);

    % Convert back to a polytope
    poly_ = mptPolytope(pZ);

%     % Visualize the result
%     hold on
%     plot(poly,[1,2],'r','Filled',true,'EdgeColor','none');
%     plot(poly_,[1,2],'b');

    % Calculate the vertices
    V_ = vertices(poly_);

    % check if the vertices are identical to the original points
    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V_',1e-12,'ByRows',true) 
          error('test_polyZonotope_mptPolytope: random test 2D (polyZonotope -> mptPolytope) failed!');
       end
    end   
end


% TEST 2-dimensional (Zonotope-Point-Case)

for j = 4:6

    % Generate random zonotope
    zono = zonotope(rand(2,j) - 0.5*ones(2,j));

    % Create point located outside the zonotope
    inter = interval(zono);
    inter = interval(zonotope([mid(inter),1.1*diag(rad(inter))]));

    ind = floor(rand() * 3.9) + 1;
    V = vertices(inter);
    p = V(:,ind);

    % Determine all vertices 
    V = vertices(zono);
    V = [V,p];

    % Construct a polytope by computation of the convex hull
    % between the zonotope and the point
    pZ1 = polyZonotope(zono);
    pZ2 = polyZonotope(p,[],[],[]);

    pZ = enclose(pZ1,pZ2);

    % Convert back to a polytope
    poly_ = mptPolytope(pZ);

%     % Visualize the result
%     hold on
%     plot(poly_,[1,2],'b');     
%     plot(V(1,:),V(2,:),'.k','MarkerSize',20);

    % Calculate the vertices
    V_ = vertices(poly_);

    % check if the vertices are identical to the original points
    for i = 1:size(V_,2)
       if ~ismembertol(V_(:,i)',V',1e-12,'ByRows',true) 
          error('test_polyZonotope_mptPolytope: random test 2D (zonotope-point-case) failed!');
       end
    end   
end



% TEST 2-dimensional (Polytope Enclosure)

for j = 1:5

    % create random polynomial zonotope
    c = (rand(2,1)-0.5*ones(2,1)) * rand()*10;
    G = rand(2,7)-0.5*ones(2,7);
    ind = datasample(1:7,4,'Replace',false);
    G(:,ind) = G(:,ind)./10;
    Grest = rand(2,1)-0.5*ones(2,1);
    expMat = [eye(2), round(rand(2,5)*2)];
    pZ = polyZonotope(c,G,Grest,expMat);

    % Over-approximate by a polytope
    poly = mptPolytope(pZ);
    
    % Extract polytope parameter
    P = get(poly,'P');
    C = P.A;
    d = P.b;
    
    % construct set of random points located inside the polynomial zonotope
    N = 10000;
    
    points = randPoint(pZ,N);
    pointsExt = randPoint(pZ,'all','extreme');
    points = [pointsExt,points];

%     % Visualize the result
%     hold on
%     plot(poly,[1,2],'b');     
%     plot(points(1,:),points(2,:),'.k');

    % check if all points are located inside the over-approximating
    % polytope
    temp = C * points - d * ones(1,size(points,2));
    
    if any(temp > 1e-12)
       error('test_polyZonotope_mptPolytope: random test 2D (polytope enclosure) failed!');
    end
end



% TEST 4-dimensional (Polynomial Zonotope -> Polytope)

for j = 10:11
    
    % Generate random polytope vertices
    points = rand(4,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind(1:min(j,length(ind))));

    % Construct a mptPolytope object from the vertices
    poly = mptPolytope(V');

    % Convert to polynomial zonotope
    pZ = polyZonotope(poly);
    
    % Convert back to a polytope
    poly_ = mptPolytope(pZ);
    
    % Calculate the vertices
    V_ = vertices(poly_);
    
    % check if the vertices are identical to the original points
    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V_',1e-12,'ByRows',true) 
          error('test_polyZonotope_mptPolytope: random test 4D (polyZonotope -> mptPolytope) failed!');
       end
    end   
end


res = 1;

%------------- END OF CODE --------------