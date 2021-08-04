function res = testLongDuration_polyZonotope_supportFunc
% testLongDuration_polyZonotope_supportFunc - unit test function for
%    calculating bounds of the polynomial zonotope along a specific
%    direction 
%
% Syntax:  
%    res = testLongDuration_polyZonotope_supportFunc
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
% Written:      30-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% RANDOM TESTS

% TEST 2-dimensional
methods = {'interval','bnb','globOpt'};

for j = 1:length(methods)
    for i = 1:3

        % create random zonotope
        c = rand(2,1)-0.5*ones(2,1);
        G = rand(2,6)-0.5*ones(2,6);
        ind = datasample(1:6,4,'Replace',false);
        G(:,ind) = G(:,ind)./10;
        Grest = rand(2,1)-0.5*ones(2,1);
        expMat = [1, round(rand(1,5)*5)];
        pZ = polyZonotope(c,G,Grest,expMat);

        % choose a random direction
        dir = rand(2,1)-0.5*ones(2,1);

        % calculate the bounds for this direction
        inter = supportFunc(pZ,dir,'range',methods{j});

        % determine random point and extreme points inside the original 
        % polynomial zonotope
        N = 10000;
        points = randPoint(pZ,N);
        pointsExt = randPoint(pZ,'all','extreme');

        points = [pointsExt,points];

        % check if the all points from the polynomial zonotope are located
        % inside the calculated range
        points_ = dir' * points;
        suc = infimum(inter) - min(points_) < 1e-10 && max(points_) - supremum(inter) < 1e-10;

        if ~suc
           error('testLongDuration_polyZonotope_supportFunc: random test 2D failed!'); 
        end
    end
end


% TEST 4-dimensional

for j = 1:length(methods)
    for i = 1:3

        % create random zonotope
        c = rand(4,1)-0.5*ones(4,1);
        G = rand(4,6)-0.5*ones(4,6);
        ind = datasample(1:6,4,'Replace',false);
        G(:,ind) = G(:,ind)./10;
        Grest = rand(4,2)-0.5*ones(4,2);
        expMat = [1, round(rand(1,5)*5)];
        pZ = polyZonotope(c,G,Grest,expMat);

        % choose a random direction
        dir = rand(4,1)-0.5*ones(4,1);

        % calculate the bounds for this direction
        inter = supportFunc(pZ,dir,'range',methods{j});

        % determine random point and extreme points inside the original 
        % polynomial zonotope
        N = 10000;
        points = randPoint(pZ,N);
        pointsExt = randPoint(pZ,'all','extreme');

        points = [pointsExt,points];

        % check if the all points from the polynomial zonotope are located
        % inside the calculated range
        points_ = dir' * points;
        suc = infimum(inter) - min(points_) < 1e-10 && max(points_) - supremum(inter) < 1e-10;

        if ~suc
           error('testLongDuration_polyZonotope_supportFunc: random test 4D failed!'); 
        end
    end
end

res = true;

%------------- END OF CODE --------------