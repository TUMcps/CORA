function res = testLong_polyZonotope_polytope
% testLong_polyZonotope_polytope - unit test function for the
%    conversion of a oolytope to a polynomial zonotope
%
% Syntax:
%    res = testLong_polyZonotope_polytope
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

% Authors:       Niklas Kochdumper
% Written:       30-October-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% TEST 2-dimensional (Polytope -> Polynomial Zonotope)

for j = 10:15
    
    % Generate random polytope vertices
    points = rand(2,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind(1:min(j,length(ind))));

    % Construct a polytope object from the vertices
    P = polytope(V);

    % Convert to polynomial zonotope
    pZres = polyZonotope(P);
    
    % check if the all polytope vertices are located inside the polynomial
    % zonotope
    suc = containsPointSet(pZres,V,[],30);
    
    if ~suc
        throw(CORAerror('CORA:testFailed'));
    end   
end


% TEST 4-dimensional (Polytope -> Polynomial Zonotope)

for j = 10:12
    
    % Generate random polytope vertices
    points = rand(4,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind(1:min(j,length(ind))));

    % Construct a polytope object from the vertices
    P = polytope(V);

    % Convert to polynomial zonotope
    pZres = polyZonotope(P);
    
    % check if the all polytope vertices are located inside the polynomial
    % zonotope
    suc = containsPointSet(pZres,V,[],30);
    
    if ~suc
        throw(CORAerror('CORA:testFailed'));
    end   
end


% TEST 2-dimensional (Polynomial Zonotope -> Polytope)

for j = 10:15

    % Generate random polytope vertices
    points = rand(2,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind(1:min(j,length(ind))));

    % Construct a polytope object from the vertices
    P = polytope(V);

    % Convert to polynomial zonotope
    pZ = polyZonotope(P);

    % Convert back to a polytope
    P_ = polytope(pZ);

%     % Visualize the result
%     hold on
%     plot(P,[1,2],'FaceColor','r');
%     plot(P_,[1,2],'b');

    % Calculate the vertices
    V_ = vertices(P_);

    % check if the vertices are identical to the original points
    if ~compareMatrices(V,V_,1e-12)
        throw(CORAerror('CORA:testFailed'));
    end   
end


% TEST 2-dimensional (Zonotope-Point-Case)

for j = 4:6

    % Generate random zonotope
    Z = zonotope(rand(2,j) - 0.5*ones(2,j));

    % Create point located outside the zonotope
    I = interval(Z);
    I = interval(zonotope([center(I),1.1*diag(rad(I))]));

    ind = floor(rand() * 3.9) + 1;
    V = vertices(I);
    p = V(:,ind);

    % Determine all vertices 
    V = vertices(Z);
    V = [V,p];

    % Construct a polytope by computation of the convex hull
    % between the zonotope and the point
    pZ1 = polyZonotope(Z);
    pZ2 = polyZonotope(p,[],[],[]);

    pZ = enclose(pZ1,pZ2);

    % Convert back to a polytope
    P_ = polytope(pZ);

    % Calculate the vertices
    V_ = vertices(P_);

%     % Visualize the result
%     hold on
%     plot(P_,[1,2],'b');     
%     plot(V(1,:),V(2,:),'.k','MarkerSize',20);
%     plot(V_(1,:),V_(2,:),'.r','MarkerSize',20);

    % check if the vertices are identical to the original points (only has
    % to be a subset)
    if ~compareMatrices(V_,V,1e-12,'subset')
        throw(CORAerror('CORA:testFailed'));
    end   
end


% TEST 2-dimensional (Polytope Enclosure)

for j = 1:5

    % create random polynomial zonotope
    c = (rand(2,1)-0.5*ones(2,1)) * rand()*10;
    G = rand(2,7)-0.5*ones(2,7);
    ind = datasample(1:7,4,'Replace',false);
    G(:,ind) = G(:,ind)./10;
    GI = rand(2,1)-0.5*ones(2,1);
    E = [eye(2), round(rand(2,5)*2)];
    pZ = polyZonotope(c,G,GI,E);

    % Over-approximate by a polytope
    P = polytope(pZ);
    
    % Extract polytope parameter
    C = P.A;
    d = P.b;
    
    % construct set of random points located inside the polynomial zonotope
    N = 10000;
    
    points = randPoint(pZ,N);
    pointsExt = randPoint(pZ,'all','extreme');
    points = [pointsExt,points];

%     % Visualize the result
%     hold on
%     plot(P,[1,2],'b');     
%     plot(points(1,:),points(2,:),'.k');

    % check if all points are located inside the over-approximating
    % polytope
    temp = C * points - d * ones(1,size(points,2));
    
    if any(temp > 1e-12)
        throw(CORAerror('CORA:testFailed'));
    end
end


% TEST 4-dimensional (Polynomial Zonotope -> Polytope)

for j = 10:11
    
    % Generate random polytope vertices
    points = rand(4,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind(1:min(j,length(ind))));

    % Construct a polytope object from the vertices
    P = polytope(V);

    % Convert to polynomial zonotope
    pZ = polyZonotope(P);
    
    % Convert back to a polytope
    P_ = polytope(pZ);
    
    % Calculate the vertices
    V_ = vertices(P_);
    
    % check if the vertices are identical to the original points
    if ~compareMatrices(V,V_,1e-12)
        throw(CORAerror('CORA:testFailed'));
    end   
end

% ------------------------------ END OF CODE ------------------------------
