function res = test_polyZonotope_mtimes
% test_polyZonotope_mtimes - unit test function for multiplication between
%                            an interval matrix and a zonotope 
%
% Syntax:  
%    res = test_polyZonotope_mtimes
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

%% ANALYTICAL TESTS

% TEST 1

% create polynomial zonotope
c = [1;2];
G = [1 2 1 -3; 1 -1 2 -1];
expMat = [1 0 0 2; 0 1 2 1];
Grest = [];
pZ = polyZonotope(c,G,Grest,expMat);

% create interval matrix
matrix = interval([1 -0.5; -1 0], [3 0.5; 3 2]);

% multiply interval matrix with polynomial zonotope
pZres = matrix * pZ;

% define ground truth
c = [2; 3];
G =  [2 4 2 -6; 2 1 3 -4];
Grest = [11.5 0; 0 23];

% check for correctness
if any(c-pZres.c) || any(any(G-pZres.G)) || any(any(Grest-pZres.Grest))
    error('test_polyZonotope_mtimes: analytical test 1 failed!');
end





%% RANDOM TESTS

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

    % create random interval matrix
    mat1 = rand(2) - 0.5*ones(2);
    mat2 = rand(2) - 0.5*ones(2);
    matrix = interval(min(mat1,mat2),max(mat1,mat2));
    
    % multiply matrix with polynomial zonotope
    pZres = matrix * pZ;

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = pointSet(pZ,N);
    pointsExt = pointSetExtreme(pZ);

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
       error('test_polyZonotope_mtimes: random test 2D failed!'); 
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

    % create random interval matrix
    mat1 = rand(4) - 0.5*ones(4);
    mat2 = rand(4) - 0.5*ones(4);
    matrix = interval(min(mat1,mat2),max(mat1,mat2));
    
    % multiply matrix with polynomial zonotope
    pZres = matrix * pZ;

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = pointSet(pZ,N);
    pointsExt = pointSetExtreme(pZ);

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
       error('test_polyZonotope_mtimes: random test 4D failed!'); 
    end
end

res = true;

%------------- END OF CODE --------------