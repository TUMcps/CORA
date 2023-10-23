function res = testLong_polyZonotope_supportFunc
% testLong_polyZonotope_supportFunc - unit test function for
%    calculating bounds of the polynomial zonotope along a specific
%    direction 
%
% Syntax:
%    res = testLong_polyZonotope_supportFunc
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
% Written:       30-July-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% TEST 2-dimensional
methods = {'interval','bnb','globOpt'};

for j = 1:length(methods)
    for i = 1:3

        % create random zonotope
        c = rand(2,1)-0.5*ones(2,1);
        G = rand(2,6)-0.5*ones(2,6);
        ind = datasample(1:6,4,'Replace',false);
        G(:,ind) = G(:,ind)./10;
        GI = rand(2,1)-0.5*ones(2,1);
        E = [1, round(rand(1,5)*5)];
        pZ = polyZonotope(c,G,GI,E);

        % choose a random direction
        dir = rand(2,1)-0.5*ones(2,1);

        % calculate the bounds for this direction
        I = supportFunc(pZ,dir,'range',methods{j});

        % determine random point and extreme points inside the original 
        % polynomial zonotope
        N = 10000;
        points = randPoint(pZ,N);
        pointsExt = randPoint(pZ,'all','extreme');

        points = [pointsExt,points];

        % check if the all points from the polynomial zonotope are located
        % inside the calculated range
        points_ = dir' * points;
        check1 = infimum(I) < min(points_) || withinTol(infimum(I) - min(points_),0,1e-10);
        check2 = max(points_) < supremum(I) || withinTol(max(points_) - supremum(I),0,1e-10);

        if ~check1 || ~check2
            res = false; return
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
        GI = rand(4,2)-0.5*ones(4,2);
        E = [1, round(rand(1,5)*5)];
        pZ = polyZonotope(c,G,GI,E);

        % choose a random direction
        dir = rand(4,1)-0.5*ones(4,1);

        % calculate the bounds for this direction
        I = supportFunc(pZ,dir,'range',methods{j});

        % determine random point and extreme points inside the original 
        % polynomial zonotope
        N = 10000;
        points = randPoint(pZ,N);
        pointsExt = randPoint(pZ,'all','extreme');

        points = [pointsExt,points];

        % check if the all points from the polynomial zonotope are located
        % inside the calculated range
        points_ = dir' * points;
        if ~(infimum(I) - min(points_) < 1e-10 && max(points_) - supremum(I) < 1e-10)
            res = false; return
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
