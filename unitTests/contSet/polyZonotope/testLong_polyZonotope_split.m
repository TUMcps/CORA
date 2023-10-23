function res = testLong_polyZonotope_split
% testLong_polyZonotope_split - unit test function for the
%    splitting of a polynomial zonotope object
%
% Syntax:
%    res = testLong_polyZonotope_split
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
    
    points1 = randPoint(pZsplit{1},N);
    pointsExt1 = randPoint(pZsplit{1},'all','extreme');
    
    points2 = randPoint(pZsplit{2},N);
    pointsExt2 = randPoint(pZsplit{2},'all','extreme');

    points = [pointsExt1,points1,pointsExt2,points2];
    
    % check if the all transformed random points are located inside the
    % resulting polynomial zonotope object
    if ~containsPointSet(pZ,points,[],30)
        res = false; return
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
    
    points1 = randPoint(pZsplit{1},N);
    pointsExt1 = randPoint(pZsplit{1},'all','extreme');
    
    points2 = randPoint(pZsplit{2},N);
    pointsExt2 = randPoint(pZsplit{2},'all','extreme');

    points = [pointsExt1,points1,pointsExt2,points2];
    
    % check if the all transformed random points are located inside the
    % resulting polynomial zonotope object
    if ~containsPointSet(pZ,points)
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
