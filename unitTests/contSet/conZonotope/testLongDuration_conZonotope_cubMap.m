function res = testLongDuration_conZonotope_cubMap
% testLongDuration_conZonotope_cubMap - unit test function for cubic
%    multiplication of constrained zonotopes
%
% Syntax:  
%    res = testLongDuration_conZonotope_cubMap
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
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;


%% RANDOM TESTS

% TEST 1: Mixed Multiplication

for i = 1:3

    % create three random zonotopes
    cZ1 = conZonotope.generateRandom(2,[],4);
    cZ2 = conZonotope.generateRandom(2,[],3);
    cZ3 = conZonotope.generateRandom(2,[],3);
    
    % create a random tensor
    T{1,1} = rand(2) - 0.5*ones(2);
    T{1,2} = rand(2) - 0.5*ones(2);
    T{2,1} = rand(2) - 0.5*ones(2);
    T{2,2} = rand(2) - 0.5*ones(2);
    
    % obtain result
    cZres = cubMap(cZ1,cZ2,cZ3,T);

    % draw random points inside the zontopes 
    N = 5;
    
    points1 = zeros(2,N);
    points2 = zeros(2,N);
    points3 = zeros(2,N);
    
    for j = 1:N
       points1(:,j) = randPoint(cZ1);
       points2(:,j) = randPoint(cZ2);
       points3(:,j) = randPoint(cZ3);
    end
    
    points1 = [points1,vertices(cZ1)];
    points2 = [points2,vertices(cZ2)];
    points3 = [points3,vertices(cZ3)];

    % calculate the cubic map for all possible point combinations
    pointsRes = zeros(2,size(points1,2)*size(points2,2)*size(points3,2));
    counter = 1;
    
    for j = 1:size(points1,2)
        for k = 1:size(points2,2)
            for h = 1:size(points3,2)
                pointsRes(:,counter) = cubMulPoint(points1(:,j),points2(:,k),points3(:,h),T);
                counter = counter + 1;
            end
        end
    end
    
%     % plot the result
%     plot(cZres,[1,2],'r','Template',50);
%     hold on
%     plot(pointsRes(1,:),pointsRes(2,:),'.k');
    
    % check if all points are located inside the calculated zonotope
    for j = 1:size(pointsRes,2)
        
       p = pointsRes(:,j);
        
       if ~in(cZres,p)
          file_name = strcat('testLongDuration_conZonotope_cubMap_1_', ...
                             datestr(now,'mm-dd-yyyy_HH-MM'));
                  
          file_path = fullfile(coraroot(), 'unitTests', 'failedTests', file_name);
          save(file_path, 'cZ1', 'cZ2', 'cZ3', 'cZres');
          
          error('conZonotope/cubMap: random test (mixed mul.) failed!'); 
       end
    end
end


% TEST 2: Cubic Multiplication

for i = 1:3

    % create random constrained zonotope
    cZ = conZonotope.generateRandom(2,[],4);

    % create a random tensor
    T{1,1} = rand(2) - 0.5*ones(2);
    T{1,2} = rand(2) - 0.5*ones(2);
    T{2,1} = rand(2) - 0.5*ones(2);
    T{2,2} = rand(2) - 0.5*ones(2);
    
    % obtain result
    cZres = cubMap(cZ,T);

    % draw random points inside the conZontope
    N = 100;
    
    points = zeros(2,N);
    
    for j = 1:N
       points(:,j) = randPoint(cZ);
    end
    
    points = [points, vertices(cZ)];

    % calculate the cubic map for all possible point combinations
    pointsRes = zeros(size(points));
    
    for j = 1:size(points,2)
        p = points(:,j);
        pointsRes(:,j) = cubMulPoint(p,p,p,T);
    end
    
%     % plot the result
%     plot(cZres,[1,2],'r','Template',50);
%     hold on
%     plot(pointsRes(1,:),pointsRes(2,:),'.k');
    
    % check if all points are located inside the calculated zonotope
    for j = 1:size(pointsRes,2)
        
       p = pointsRes(:,j);
        
       if ~in(cZres,p)
          file_name = strcat('testLongDuration_conZonotope_cubMap_2_', ...
                             datestr(now,'mm-dd-yyyy_HH-MM'));
                  
          file_path = fullfile(coraroot(), 'unitTests', 'failedTests', file_name);
          save(file_path, 'cZ', 'cZres')
          
          error('conZonotope/cubMap: random test failed!'); 
       end
    end
end

res = true;

end


% Auxiliary Functions -----------------------------------------------------

function p = cubMulPoint(x1,x2,x3,T)

    p = zeros(size(x1));
    
    % loop over all dimensions
    for i = 1:length(p)
       
         % loop over all quadratic matrices for this dimension
         for j = 1:size(T,2)
            p(i) = p(i) + (x1' * T{i,j} * x2) * x3(j);
         end
    end
end

%------------- END OF CODE --------------