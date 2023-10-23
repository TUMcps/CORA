function res = testLong_polyZonotope_mtimes
% testLong_polyZonotope_mtimes - unit test function for
%    multiplication between an interval matrix and a zonotope 
%
% Syntax:
%    res = testLong_polyZonotope_mtimes
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

    % create random interval matrix
    mat1 = rand(2) - 0.5*ones(2);
    mat2 = rand(2) - 0.5*ones(2);
    matrix = interval(min(mat1,mat2),max(mat1,mat2));
    
    % multiply matrix with polynomial zonotope
    pZres = matrix * pZ;

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = randPoint(pZ,N);
    pointsExt = randPoint(pZ,'all','extreme');

    points = [pointsExt,points];
    
    % determine random sample matrices from the interval matrix
    N = size(points,2);
    mat = cell(N,1);
    rMat = rad(matrix);
    cMat = center(matrix);
    
    for j = 1:N
        temp = rand(2)*2 - ones(2);
        mat{j} = cMat + temp.*rMat;
    end
    
    % multiply the random matrices with the random points from the original
    % polynomial zonotope
    points_ = zeros(2,N);
    
    for j = 1:N
       points_(:,j) = mat{j} * points(:,j);
    end

    % check if the all transformed random points are located inside the
    % resulting polynomial zonotope object
    suc = containsPointSet(pZres,points_,[],30);
    
    if ~suc
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

    % create random interval matrix
    mat1 = rand(4) - 0.5*ones(4);
    mat2 = rand(4) - 0.5*ones(4);
    matrix = interval(min(mat1,mat2),max(mat1,mat2));
    
    % multiply matrix with polynomial zonotope
    pZres = matrix * pZ;

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = randPoint(pZ,N);
    pointsExt = randPoint(pZ,'all','extreme');

    points = [pointsExt,points];
    
    % determine random sample matrices from the interval matrix
    N = size(points,2);
    mat = cell(N,1);
    rMat = rad(matrix);
    cMat = center(matrix);
    
    for j = 1:N
        temp = rand(4)*2 - ones(4);
        mat{j} = cMat + temp.*rMat;
    end
    
    % multiply the random matrices with the random points from the original
    % polynomial zonotope
    points_ = zeros(4,N);
    
    for j = 1:N
       points_(:,j) = mat{j} * points(:,j);
    end

    % check if the all transformed random points are located inside the
    % resulting polynomial zonotope object
    suc = containsPointSet(pZres,points_);
    
    if ~suc
       throw(CORAerror('CORA:testFailed'));
    end
end

% ------------------------------ END OF CODE ------------------------------
