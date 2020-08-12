function res = test_polyZonotope_zonotope
% test_polyZonotope_zonotope - unit test function for zonotope
%                              over-approximation of a polynomial zonotope
%
% Syntax:  
%    res = test_polyZonotope_zonotope
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

% reduce the polynomial zonotope
zono = zonotope(pZ);
c_ = zono.Z(:,1);
G_ = zono.Z(:,2:end);

% define ground truth
c = [1.5; 3];
G =  [1 2 0.5 -3; 1 -1 1 -1];

% check for correctness
if any(c-c_)
    error('test_polyZonotope_zonotope: analytical test 1 failed!');
end
for i = 1:size(G,2)
   if ~ismember(G(:,i)',G_)
       error('test_polyZonotope_zonotope: analytical test 1 failed!');
   end
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

    % calculate zonotope over-approximation
    zono = zonotope(pZ);
    zono = halfspace(zono);

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = pointSet(pZ,N);
    pointsExt = pointSetExtreme(pZ);

    points = [pointsExt,points];

    % check if the all points from the original polynomial zonotope are
    % enclosed by the reduced polynomial zonotope
    suc = 1;
    for j = 1:size(points,2)
       if ~containsPoint(zono,points(:,j),1e-10)
           suc = 0;
           break;
       end
    end
    
    if ~suc
       error('test_polyZonotope_zonotope: random test 2D failed!'); 
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

    % calculate zonotope over-approximation
    zono = zonotope(pZ);
    zono = halfspace(zono);

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = pointSet(pZ,N);
    pointsExt = pointSetExtreme(pZ);

    points = [pointsExt,points];

    % check if the all points from the original polynomial zonotope are
    % enclosed by the reduced polynomial zonotope
    suc = 1;
    for j = 1:size(points,2)
       if ~containsPoint(zono,points(:,j),1e-10)
           suc = 0;
           break;
       end
    end
    
    if ~suc
       error('test_polyZonotope_zonotope: random test 4D failed!'); 
    end
end

res = true;

%------------- END OF CODE --------------