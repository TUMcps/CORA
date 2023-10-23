function res = testLong_polyZonotope_plus
% testLong_polyZonotope_plus - unit test function for the addition
%    of two polynomial zonotope objects
%
% Syntax:
%    res = testLong_polyZonotope_plus
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
% Written:       26-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% TEST linear zonotopes

% express linear zonotopes as polynomial zonotopes and compare the result
% of addition to the result for addition of linear zonotopes

% loop over different dimensions
for i = 2:10

    % create random linear zonotopes
    G1 = rand(i,5);
    c1 = rand(i,1);

    G2 = rand(i,7);
    c2 = rand(i,1);

    Z1 = zonotope(c1,G1);
    Z2 = zonotope(c2,G2);

    % addition of the linear zonotopes
    Z_ = Z1 + Z2;

    % express lienar zonotopes as polynomial zonotopes
    E1 = eye(size(G1,2));
    E2 = [zeros(size(G1,2),size(G2,2)); eye(size(G2,2))];

    pZ1 = polyZonotope(c1,G1,[],E1);
    pZ2 = polyZonotope(c2,G2,[],E2);

    % addition of the polynomial zonotopes
    pZres = pZ1 + pZ2;
    
    % check for correctness of the result
    c_ = Z_.c;
    G_ = Z_.G;
    
    if ~compareMatrices(pZres.c,c_) || ~compareMatrices(pZres.G,G_)
        throw(CORAerror('CORA:testFailed'));
    end
    
end


% TEST 2-dimensional

for i = 1:5
    
    % create random zonotope
    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,7)-0.5*ones(2,7);
    ind = datasample(1:7,4,'Replace',false);
    G(:,ind) = G(:,ind)./10;
    GI = rand(2,1)-0.5*ones(2,1);
    E = [eye(2), round(rand(2,5)*5)];
    pZ = polyZonotope(c,G,GI,E);

    % create random transformation matrix
    matrix = rand(2) - 0.5*ones(2);
    
    % add the polynomial zonotope and the transformed polynomial zonotope
    pZres = exactPlus(pZ,matrix * pZ);

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = randPoint(pZ,N);
    pointsExt = randPoint(pZ,'all','extreme');

    points = [pointsExt,points];
    
    % transform the random points with the transformation matrix
    N = size(points,2);
    
    for j = 1:N
       points(:,j) = points(:,j) + matrix * points(:,j); 
    end
    
    % check if the all transformed random points are located inside the
    % resulting polynomial zonotope object
    if ~containsPointSet(pZres,points,[],30)
        throw(CORAerror('CORA:testFailed'));
    end
end


% TEST 4-dimensional

for i = 1:5
    
    % create random zonotope
    c = rand(4,1)-0.5*ones(4,1);
    G = rand(4,6)-0.5*ones(4,6);
    ind = datasample(1:6,4,'Replace',false);
    G(:,ind) = G(:,ind)./10;
    GI = rand(4,2)-0.5*ones(4,2);
    E = [eye(4), round(rand(4,2)*5)];
    pZ = polyZonotope(c,G,GI,E);

    % create random transformation matrix
    matrix = rand(4) - 0.5*ones(4);
    
    % add the polynomial zonotope and the transformed polynomial zonotope
    pZres = exactPlus(pZ, matrix * pZ);

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = randPoint(pZ,N);
    pointsExt = randPoint(pZ,'all','extreme');

    points = [pointsExt,points];
    
    % transform the random points with the transformation matrix
    N = size(points,2);
    
    for j = 1:N
        points(:,j) = points(:,j) + matrix * points(:,j); 
    end
    
    % check if the all transformed random points are located inside the
    % resulting polynomial zonotope object
    if ~containsPointSet(pZres,points)
       throw(CORAerror('CORA:testFailed'));
    end
end

% ------------------------------ END OF CODE ------------------------------
