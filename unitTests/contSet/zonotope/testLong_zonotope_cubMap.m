function res = testLong_zonotope_cubMap
% testLong_zonotope_cubMap - unit test function for cubic multiplication of 
%                        zonotopes
%
% Syntax:
%    res = testLong_zonotope_cubMap
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
% Written:       16-August-2018
% Last update:   01-May-2020 (MW, cubicMultiplication -> cubMap)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TEST 1: Mixed Multiplication

for i = 1:10

    % create three random zonotopes
    Z1 = zonotope(rand(2,3)-0.5*ones(2,3));
    Z2 = zonotope(rand(2,4)-0.5*ones(2,4));
    Z3 = zonotope(rand(2,5)-0.5*ones(2,5));
    
    % create a random tensor
    T{1,1} = rand(2) - 0.5*ones(2);
    T{1,2} = rand(2) - 0.5*ones(2);
    T{2,1} = rand(2) - 0.5*ones(2);
    T{2,2} = rand(2) - 0.5*ones(2);
    
    % obtain result
    Zres = cubMap(Z1,Z2,Z3,T);

    % draw random points inside the zonotopes 
    N = 20;
    
    points1 = zeros(2,N);
    points2 = zeros(2,N);
    points3 = zeros(2,N);
    
    for j = 1:N
       points1(:,j) = randPoint(Z1);
       points2(:,j) = randPoint(Z2);
       points3(:,j) = randPoint(Z3);
    end
    
    points1 = [points1,vertices(Z1)];
    points2 = [points2,vertices(Z2)];
    points3 = [points3,vertices(Z3)];

    % calculate the cubic map for all possible point combinations
    pointsRes = zeros(2,size(points1,2)*size(points2,2)*size(points3,2));
    counter = 1;
    
    for j = 1:size(points1,2)
        for k = 1:size(points2,2)
            for h = 1:size(points3,2)
                pointsRes(:,counter) = cubMapPoint(points1(:,j),points2(:,k),points3(:,h),T);
                counter = counter + 1;
            end
        end
    end
    
    % convert zonotope to halfspace representation
    Zres = halfspace(Zres);
    C = Zres.halfspace.H;
    d = Zres.halfspace.K;
    
%     % plot the result
%     plot(Zres,[1,2],'r');
%     hold on
%     plot(pointsRes(1,:),pointsRes(2,:),'.k');
    
    % check if all points are located inside the calculated zonotope
    for j = 1:size(pointsRes,2)
        
       p = pointsRes(:,j);
       if ~all(C*p - d < 0 | withinTol(C*p - d,0,1e-12))
          throw(CORAerror('CORA:testFailed'));
       end
    end
end


% TEST 2: Cubic Multiplication

for i = 1:10

    % create random zonotope
    Z = zonotope(rand(2,4)-0.5*ones(2,4));
    
    % create a random tensor
    T{1,1} = rand(2) - 0.5*ones(2);
    T{1,2} = rand(2) - 0.5*ones(2);
    T{2,1} = rand(2) - 0.5*ones(2);
    T{2,2} = rand(2) - 0.5*ones(2);
    
    % obtain result
    Zres = cubMap(Z,T);

    % draw random points inside the zonotopes 
    N = 1000;
    
    points = zeros(2,N);
    
    for j = 1:N
       points(:,j) = randPoint(Z);
    end
    
    points = [points, vertices(Z)];

    % calculate the cubic map for all possible point combinations
    pointsRes = zeros(size(points));
    
    for j = 1:size(points,2)
        p = points(:,j);
        pointsRes(:,j) = cubMapPoint(p,p,p,T);
    end
    
    % convert zonotope to halfspace representation
    Zres = halfspace(Zres);
    C = Zres.halfspace.H;
    d = Zres.halfspace.K;
    
%     % plot the result
%     plot(Zres,[1,2],'r');
%     hold on
%     plot(pointsRes(1,:),pointsRes(2,:),'.k');
    
    % check if all points are located inside the calculated zonotope
    for j = 1:size(pointsRes,2)
        
       p = pointsRes(:,j);
        
       if ~all(C*p - d < 0 | withinTol(C*p - d,0,1e-12))
          throw(CORAerror('CORA:testFailed'));
       end
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
