function res = test_polyZonotope_reduce
% test_polyZonotope_reduce - unit test function of order reduction
%
% Syntax:  
%    res = test_polyZonotope_reduce
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
% Written:      29-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% TEST 1

% create polynomial zonotope
c = [0;0];
G = [0 2 1;3 0 1];
expMat = [1 0 1;0 1 1];
Grest = [0;0];
pZ = polyZonotope(c,G,Grest,expMat);

% reduce the polynomial zonotope
pZred = reduce(pZ,'girard',1);

% define ground truth
Grest =  [3 0; 0 4];
c = [0;0];

% check for correctness
if any(c-pZred.c) || any(any(Grest-pZred.Grest)) || ~isempty(pZred.G) || ...
   ~isempty(pZred.expMat) || ~isempty(pZred.id)
    error('test_polyZonotope_reduce: analytical test 1 failed!');
end



% TEST 2

% create polynomial zonotope
c = [0;0];
G = [0 3 1 1; 2 0 1 -3];
expMat = [1 0 1 2;0 1 1 0];
Grest = [-1 -4; -2 -1];
pZ = polyZonotope(c,G,Grest,expMat);

% reduce the polynomial zonotope
pZred = reduce(pZ,'girard',2);

% define ground truth
c = [0.5;-1.5];
G = [3;0];
Grest = [-4 2.5 0; -1 0 6.5];
expMat = 1;

% check for correctness
if any(c-pZred.c) || any(any(G-pZred.G)) || any(any(Grest-pZred.Grest)) || ...
   any(any(expMat-pZred.expMat))
    error('test_polyZonotope_reduce: analytical test 2 failed!');
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

    % reduce the polynomial zonotope
    pZred = reduce(pZ,'girard',2);

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = pointSet(pZ,N);
    pointsExt = pointSetExtreme(pZ);

    points = [pointsExt,points];

    % check if the all points from the original polynomial zonotope are
    % enclosed by the reduced polynomial zonotope
    suc = containsPointSet(pZred,points,[],30);
    
    if ~suc
       error('test_polyZonotope_reduce: random test 2D failed!'); 
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

    % reduce the polynomial zonotope
    pZred = reduce(pZ,'girard',2);

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = pointSet(pZ,N);
    pointsExt = pointSetExtreme(pZ);

    points = [pointsExt,points];

    % check if the all points from the original polynomial zonotope are
    % enclosed by the reduced polynomial zonotope
    suc = containsPointSet(pZred,points);
    
    if ~suc
       error('test_polyZonotope_reduce: random test 4D failed!'); 
    end
end

res = true;

%------------- END OF CODE --------------