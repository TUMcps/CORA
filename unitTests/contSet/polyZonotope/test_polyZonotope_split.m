function res = test_polyZonotope_split
% test_polyZonotope_split - unit test function for the splitting of a
%                           polynomial zonotope object
%
% Syntax:  
%    res = test_polyZonotope_split
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

% create polynomial zonotopes
c = [-1;3];
G = [2 0 1; 1 2 1];
expMat = [1 0 2; 0 1 1];
Grest = [];
pZ = polyZonotope(c,G,Grest,expMat);

% plotRandPoint(pZ,[1,2],100000,'.r');

% split the polynomial zonotope at the longest generator
pZsplit = splitLongestGen(pZ);

% define ground truth
c1 = [0; 3.5];
G1 =  [1 1/4 1/2 1/4; 1/2 9/4 1/2 1/4];
expMat1 = [1 0 1 2; 0 1 1 1];

c2 = [-2; 2.5];
G2 =  [1 1/4 -1/2 1/4; 1/2 9/4 -1/2 1/4];
expMat2 = [1 0 1 2; 0 1 1 1];

% check for correctness
if any(c1-pZsplit{1}.c) || any(c2-pZsplit{2}.c)
    error('test_polyZonotope_split: analytical test 1 failed!');
end

for i = 1:size(expMat1,2)    
    
    ind = ismember(pZsplit{1}.expMat',expMat1(:,i)','rows');  
    ind_ = find(ind > 0);
    
    if isempty(ind_)
        error('test_polyZonotope_split: analytical test failed!');        
    elseif ~all(pZsplit{1}.G(:,ind_(1)) == G1(:,i))
        error('test_polyZonotope_split: analytical test failed!');
    end
end

for i = 1:size(expMat2,2)    
    
    ind = ismember(pZsplit{2}.expMat',expMat2(:,i)','rows');  
    ind_ = find(ind > 0);
    
    if isempty(ind_)
        error('test_polyZonotope_split: analytical test failed!');        
    elseif ~all(pZsplit{2}.G(:,ind_(1)) == G2(:,i))
        error('test_polyZonotope_split: analytical test failed!');
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

    % split the polynomial zonotope at the longest generator
    if i < 3
        pZsplit = splitLongestGen(pZ);
    else
        polyOrd = 2;
        pZsplit = splitLongestGen(pZ,polyOrd);
    end

    % determine random point and extreme points inside the splitted 
    % polynomial zonotopes
    N = 5000;
    
    points1 = pointSet(pZsplit{1},N);
    pointsExt1 = pointSetExtreme(pZsplit{1});
    
    points2 = pointSet(pZsplit{2},N);
    pointsExt2 = pointSetExtreme(pZsplit{2});

    points = [pointsExt1,points1,pointsExt2,points2];
    
    % check if the all transformed random points are located inside the
    % resulting polynomial zonotope object
    suc = containsPointSet(pZ,points,[],30);
    
    if ~suc
       error('test_polyZonotope_split: random test 2D failed!'); 
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

    % split the polynomial zonotope at the longest generator
    if i < 3
        pZsplit = splitLongestGen(pZ);
    else
        polyOrd = 2;
        pZsplit = splitLongestGen(pZ,polyOrd);
    end

    % determine random point and extreme points inside the splitted 
    % polynomial zonotopes
    N = 5000;
    
    points1 = pointSet(pZsplit{1},N);
    pointsExt1 = pointSetExtreme(pZsplit{1});
    
    points2 = pointSet(pZsplit{2},N);
    pointsExt2 = pointSetExtreme(pZsplit{2});

    points = [pointsExt1,points1,pointsExt2,points2];
    
    % check if the all transformed random points are located inside the
    % resulting polynomial zonotope object
    suc = containsPointSet(pZ,points);
    
    if ~suc
       error('test_polyZonotope_split: random test 4D failed!'); 
    end
end

res = true;

%------------- END OF CODE --------------