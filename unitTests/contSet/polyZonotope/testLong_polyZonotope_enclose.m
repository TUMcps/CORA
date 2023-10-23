function res = testLong_polyZonotope_enclose
% testLong_polyZonotope_enclose - unit test function for
%    convex-hull computation of polynomial zonotopes
%
% Syntax:
%    res = testLong_polyZonotope_enclose
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
% Written:       19-July-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TEST 2-dimensional

for i = 1:5
    
    % create random polynomial zonotope 1
    c = (rand(2,1)-0.5*ones(2,1)) * rand()*10;
    G = rand(2,7)-0.5*ones(2,7);
    ind = datasample(1:7,4,'Replace',false);
    G(:,ind) = G(:,ind)./10;
    GI = rand(2,1)-0.5*ones(2,1);
    E = [eye(2), round(rand(2,5)*5)];
    pZ1 = polyZonotope(c,G,GI,E);
    
    % create a random linear transformation and apply it to the polynomial
    % zonotope
    A = (rand(2)-0.5*ones(2)) * 5;
    b = (rand(2,1)-0.5*ones(2,1)) * 10;
    
    pZ2 = A*pZ1 + b;
    
    % create random zonotope and add it to the polynomial zonotope
    zono = zonotope(rand(2,3)-0.5*ones(2,3));
    
    pZ2 = pZ2 + zono;

    % compute the convex hull
    pZ = enclose(pZ1,pZ2);

    % determine random point and extreme points inside the original two 
    % polynomial zonotopes
    N = 1000;
    H = 10;
    
    points = randPoint(pZ1,N);
    pointsExt = randPoint(pZ1,'all','extreme');
    points1 = [pointsExt,points];
    
    points2 = A*points1 + b*ones(1,size(points1,2));
    
    % calculate linear combinations between all points of the two
    % polynomial zonotopes
    points = zeros(size(points1,1),N*H);
    space = linspace(0,1,H);
    counter = 1;
   
    for j = 1:N
        for h = 1:10
           points(:,counter) = space(h) * points1(:,j) + ...
                               (1-space(h)) * (points2(:,j) + randPoint(zono));
           counter = counter + 1;
        end
    end
    
%     % visualize result
%     figure; hold on;
%     plot(pZ,[1,2],'FaceColor',[.6 .6 .6],'Order',3);
%     plotRandPoint(pZ1,[1,2],100000,'.r');
%     plotRandPoint(pZ2,[1,2],100000,'.b');
%     plot(points(1,:),points(2,:),'.k');

    % check if the all points from the original polynomial zonotope are
    % enclosed by the reduced polynomial zonotope
    if ~containsPointSet(pZ,points,[],30)
        throw(CORAerror('CORA:testFailed'));
    end
end


% TEST 4-dimensional

for i = 1:5
    
    % create random polynomial zonotope 1
    c = (rand(4,1)-0.5*ones(4,1)) * rand()*10;
    G = rand(4,7)-0.5*ones(4,7);
    ind = datasample(1:7,4,'Replace',false);
    G(:,ind) = G(:,ind)./10;
    GI = rand(4,1)-0.5*ones(4,1);
    E = [eye(4), round(rand(4,3)*5)];
    pZ1 = polyZonotope(c,G,GI,E);
    
    % create a random linear transformation and apply it to the polynomial
    % zonotope
    A = (rand(4)-0.5*ones(4)) * 5;
    b = (rand(4,1)-0.5*ones(4,1)) * 10;
    
    pZ2 = A*pZ1 + b;
    
    % create random zonotope and add it to the polynomial zonotope
    zono = zonotope(rand(4,3)-0.5*ones(4,3));
    
    pZ2 = pZ2 + zono;

    % compute the convex hull
    pZ = enclose(pZ1,pZ2);

    % determine random point and extreme points inside the original two 
    % polynomial zonotopes
    N = 1000;
    H = 10;
    
    points = randPoint(pZ1,N);
    pointsExt = randPoint(pZ1,'all','extreme');
    points1 = [pointsExt,points];
    
    points2 = A*points1 + b*ones(1,size(points1,2));
    
    % calculate linear combinations between all points of the two
    % polynomial zonotopes
    points = zeros(size(points1,1),N*H);
    space = linspace(0,1,H);
    counter = 1;
   
    for j = 1:N
        for h = 1:10
           points(:,counter) = space(h) * points1(:,j) + ...
                               (1-space(h)) * (points2(:,j) + randPoint(zono));
           counter = counter + 1;
        end
    end

    % check if the all points from the original polynomial zonotope are
    % enclosed by the reduced polynomial zonotope
    if ~containsPointSet(pZ,points,[],5)
       throw(CORAerror('CORA:testFailed'));
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
