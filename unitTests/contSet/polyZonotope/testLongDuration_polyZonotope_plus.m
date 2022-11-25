function res = testLongDuration_polyZonotope_plus
% testLongDuration_polyZonotope_plus - unit test function for the addition
%    of two polynomial zonotope objects
%
% Syntax:  
%    res = test_polyZonotope_plus
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
% Written:      26-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% RANDOM TESTS


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

    zono1 = zonotope([c1,G1]);
    zono2 = zonotope([c2,G2]);

    % addition of the linear zonotopes
    zonoRes = zono1 + zono2;

    % express lienar zonotopes as polynomial zonotopes
    expMat1 = eye(size(G1,2));
    expMat2 = [zeros(size(G1,2),size(G2,2)); eye(size(G2,2))];

    pZ1 = polyZonotope(c1,G1,[],expMat1);
    pZ2 = polyZonotope(c2,G2,[],expMat2);

    % addition of the polynomial zonotopes
    pZres = pZ1 + pZ2;
    
    % check for correctness of the result
    c_ = zonoRes.Z(:,1);
    G_ = zonoRes.Z(:,2:end);
    
    if any(abs(pZres.c-c_) > 1e-15) || any(any(abs(pZres.G-G_) > 1e-15))
        error('test_polyZonotope_plus: random test linear zonotopes failed!'); 
    end
    
end



% TEST 2-dimensional

for i = 1:5
    
    % create random zonotope
    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,7)-0.5*ones(2,7);
    ind = datasample(1:7,4,'Replace',false);
    G(:,ind) = G(:,ind)./10;
    Grest = rand(2,1)-0.5*ones(2,1);
    expMat = [eye(2), round(rand(2,5)*5)];
    pZ = polyZonotope(c,G,Grest,expMat);

    % create random transformation matrix
    matrix = rand(2) - 0.5*ones(2);
    
    % add the polynomial zontope and the transformed polynomial zonotope
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
    suc = containsPointSet(pZres,points,[],30);
    
    if ~suc
       error('test_polyZonotope_plus: random test 2D failed!'); 
    end
end


% TEST 4-dimensional

for i = 1:5
    
    % create random zonotope
    c = rand(4,1)-0.5*ones(4,1);
    G = rand(4,6)-0.5*ones(4,6);
    ind = datasample(1:6,4,'Replace',false);
    G(:,ind) = G(:,ind)./10;
    Grest = rand(4,2)-0.5*ones(4,2);
    expMat = [eye(4), round(rand(4,2)*5)];
    pZ = polyZonotope(c,G,Grest,expMat);

    % create random transformation matrix
    matrix = rand(4) - 0.5*ones(4);
    
    % add the polynomial zontope and the transformed polynomial zonotope
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
    suc = containsPointSet(pZres,points);
    
    if ~suc
       error('test_polyZonotope_plus: random test 4D failed!'); 
    end
end

res = true;

%------------- END OF CODE --------------