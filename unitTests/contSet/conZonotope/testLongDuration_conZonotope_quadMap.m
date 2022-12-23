function res = testLongDuration_conZonotope_quadMap
% testLongDuration_conZonotope_quadMap - unit test function for the
%    computation of the quadratic map of a constrained zonotope
%
% Syntax:  
%    res = testLongDuration_conZonotope_quadMap
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
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      13-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;


% TEST 1: Random Test (quadratic multiplication) --------------------------

for k = 1:3
    
    % Generate random conZonotope object
    V = rand(2,3)-0.5*ones(2,3);
    P = mptPolytope(V');
    cZ = conZonotope(P);

    % generate quadratic mapping matrices
    Q{1} = rand(2)-0.5*ones(2);
    Q{2} = rand(2)-0.5*ones(2);

    % quadratic multiplication 
    cZquad = quadMap(cZ,Q);
    
    % create random points inside the original constrained zonotope
    N = 1000;
    points = zeros(2,N);
    
    for i = 1:size(points,2)
        temp = rand(3,1);
        temp = temp/sum(temp);
        points(:,i) = V(:,1)*temp(1) + V(:,2)*temp(2) + V(:,3)*temp(3);
    end
    
    points = [points,V];

    % calculate points that have to be located inside the resulting 
    % conZonotope
    points_ = zeros(size(points));

    for i = 1:size(points,2)
       p = points(:,i);
       points_(1,i) = p' * Q{1} * p;
       points_(2,i) = p' * Q{2} * p;
    end

    % convert the resulting conZonotope to a mptPolytope (to easily check if
    % a point is located inside the conZonotope)
    P = mptPolytope(cZquad);

    % extract inequality constraints
    temp = get(P,'P');
    A = temp.A;
    b = temp.b;

%     % visualize result
%     plot(points_(1,:),points_(2,:),'.k');
%     hold on
%     plot(cZquad,[1,2],'r');

    % check if all points are located inside the resulting conZonotope
    for i = 1:size(points_,2)
       if ~all(A*points_(:,i) < b | withinTol(A*points_(:,i),b))
           file_name = strcat('testLongDuration_conZonotope_quadMap_1_', ...
                              datestr(now,'mm-dd-yyyy_HH-MM'));
                  
           file_path = fullfile(CORAROOT, 'unitTests', 'failedTests', ...
                                file_name);
                           
           save(file_path, 'V','Q')
           throw(CORAerror('CORA:testFailed'));
       end
    end
end


% TEST 2: Random Test (quadratic multiplication) --------------------------

% compare the results with the implementation of the
% quadMap-function for zonotopes. If the conZonotope
% object has no constraints, the result has to be identical

for k = 2:5
   
    % create random zonotope
    Z = zonotope(rand(k,10)-0.5*ones(k,10));
    cZ = conZonotope(Z.Z);
    
    % generate random quadratic mapping matrices
    Q{1} = rand(k)-0.5*ones(k);
    Q{2} = rand(k)-0.5*ones(k);
    
    % compute quadratic map
    Zquad = quadMap(Z,Q);
    cZquad = quadMap(cZ,Q);
    
    % compare the results
    if ~compareMatrices(Zquad.Z,cZquad.Z) || ...
        ~isempty(cZquad.A) || ~isempty(cZquad.b)
        file_name = strcat('testLongDuration_conZonotope_quadMap_2_', ...
                           datestr(now,'mm-dd-yyyy_HH-MM'));
                  
        file_path = fullfile(CORAROOT, 'unitTests', 'failedTests', ...
                             file_name);
                           
        save(file_path, 'Z', 'cZ')
        throw(CORAerror('CORA:testFailed'));
    end  
end



% TEST 3: Random Test (mixed multiplication) ------------------------------

for k = 1:3
    
    % Generate first random conZonotope object
    V1 = rand(2,3)-0.5*ones(2,3);
    P = mptPolytope(V1');
    cZ1 = conZonotope(P);
    
    % Generate first random conZonotope object
    V2 = rand(2,2)-0.5*ones(2,2);
    Z = zonotope(interval(min(V2,[],2),max(V2,[],2)));
    cZ2 = conZonotope(Z.Z,[1 -1],0);
    V2 = [max(V2,[],2),min(V2,[],2)];

    % generate random quadratic mapping matrices
    Q{1} = rand(2)-0.5*ones(2);
    Q{2} = rand(2)-0.5*ones(2);

    % quadratic multiplication 
    cZquad = quadMap(cZ1,cZ2,Q);
    
    % create random points inside the original constrained zonotopes
    N = 100;
    points1 = zeros(2,N);
    points2 = zeros(2,N);
    
    for i = 1:N
        
        % first constrained zonotope
        temp = rand(3,1);
        temp = temp/sum(temp);
        points1(:,i) = V1(:,1)*temp(1) + V1(:,2)*temp(2) + V1(:,3)*temp(3);
        
        % second constrained zonotope
        temp = rand(2,1);
        temp = temp/sum(temp);
        points2(:,i) = V2(:,1)*temp(1) + V2(:,2)*temp(2);
    end
    
    points1 = [points1,V1];
    points2 = [points2,V2];

    % calculate points that have to be located inside the resulting 
    % conZonotope
    points_ = zeros(2,N^2);
    counter = 1;

    for i = 1:size(points1,2)
        for j = 1:size(points2,2)
           p1 = points1(:,i);
           p2 = points2(:,j);
           points_(1,counter) = p1' * Q{1} * p2;
           points_(2,counter) = p1' * Q{2} * p2;
           counter = counter + 1;
        end
    end

    % convert the resulting conZonotope to a mptPolytope (to easily check if
    % a point is located inside the conZonotope)
    P = mptPolytope(cZquad);

    % extract inequality constraints
    temp = get(P,'P');
    A = temp.A;
    b = temp.b;

%     % visualize result
%     plot(points_(1,:),points_(2,:),'.k');
%     hold on
%     plot(cZquad,[1,2],'r');

    % check if all points are located inside the resulting conZonotope
    for i = 1:size(points_,2)
       if ~all(A*points_(:,i) < b | withinTol(A*points_(:,i),b))
           file_name = strcat('testLongDuration_conZonotope_quadMap_3_', ...
                              datestr(now,'mm-dd-yyyy_HH-MM'));
                  
           file_path = fullfile(CORAROOT, 'unitTests', 'failedTests', ...
                                file_name);
                           
           save(file_path, 'temp')
           throw(CORAerror('CORA:testFailed'));
       end
    end
end

%------------- END OF CODE --------------