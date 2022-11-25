function res = testLongDuration_polyZonotope_cubMap
% testLongDuration_polyZonotope_cubMap - unit test function for the cubic 
%    multiplication of polynomial zonotopes with a tensor
%
% Syntax:  
%    res = testLongDuration_polyZonotope_cubMap
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
% Written:      17-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% RANDOM TESTS

% TEST zonotope (cubic multiplication)

% compare the result of the polyZonotope/cubMap function
% with the one from the zonotope/cubMap function

for dim = 2:5
   
    % create random zonotope
    zono = zonotope(rand(2,10)-0.5*ones(2,10));
    pZ = polyZonotope(zono);
    
    % create random third-order tensor
    T{1,1} = rand(2)-0.5*ones(2);
    T{1,2} = rand(2)-0.5*ones(2);
    T{2,1} = rand(2)-0.5*ones(2);
    T{2,2} = rand(2)-0.5*ones(2);
    
    % compute cubic map
    zonoRes = cubMap(zono,T);
    pZres = cubMap(pZ,T);
    
    % over approximate the resulting polynomial zonotope with a zonotope
    zonoRes_ = zonotope(pZres);
    
    % test for equality
    Tol = 1e-12;
    
    if any(size(zonoRes_) ~= size(zonoRes_)) || any(abs(center(zonoRes_) - center(zonoRes)) > Tol)
        error('testLongDuration_polyZonotope_cubMap: comparison with zonotope implementation failed!'); 
    else
        G = generators(zonoRes);
        G_ = generators(zonoRes_);
        for i = 1:size(G,2)
            if ~ismembertol(G(:,i)',G_',Tol,'ByRows',true)
                error('testLongDuration_polyZonotope_cubMap: comparison with zonotope implementation failed!'); 
            end
        end
    end
end


% TEST 2-dimensional (cubic multiplication)

for i = 1:3

    % create random zonotope
    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,8)-0.5*ones(2,8);
    Grest = rand(2,3)-0.5*ones(2,3);
    expMat = [eye(3), round(rand(3,5)*10)];
    pZ = polyZonotope(c,G,Grest,expMat);

    % create random third-order tensor
    T{1,1} = rand(2)-0.5*ones(2);
    T{1,2} = rand(2)-0.5*ones(2);
    T{2,1} = rand(2)-0.5*ones(2);
    T{2,2} = rand(2)-0.5*ones(2);

    % obtain result
    pZres = cubMap(pZ,T);

    % draw random points from zonotope and approximately compute the
    % cubic mapping
    N = 10000;
    points = zeros(2,N);

    initPoints = [randPoint(pZ,N),randPoint(pZ,'all','extreme')];

    for j = 1:size(initPoints,2)
        p = initPoints(:,j);
        points(:,j) = cubMulPoint(p,p,p,T);
    end

%     % visualize result
%     plot(pZ,[1,2],5,'FaceColor',[.5 .5 .5],'Filled',true,'Splits',30);
%     hold on
%     plot(initPoints(1,:),initPoints(2,:),'.k');
%     plot(zonotope(pZ),[1,2],'r');
% 
%     figure
%     plot(pZres,[1,2],'FaceColor',[.5 .5 .5],'Filled',true,'Order',5,'Splits',30)
%     hold on
%     plot(zonotope(pZres),[1,2],'r');
%     plot(points(1,:),points(2,:),'.k');

    % check if the polynomial zonotope obtained by computing the quadratic map
    % encloses all randomly drawn and mapped points
    suc = containsPointSet(pZres,points,[],30);
   
    if ~suc
        error('testLongDuration_polyZonotope_cubMap: random test cubic multiplication failed!');
    end  
end


% TEST zonotope (mixed cubic multiplication)

% compare the result of the polyZonotope/cubMap function
% with the one from the zonotope/cubMap function

for dim = 2:5
   
    % create random zonotopes
    zono1 = zonotope(rand(2,3)-0.5*ones(2,3));
    pZ1 = polyZonotope(zono1);
    
    zono2 = zonotope(rand(2,4)-0.5*ones(2,4));
    pZ2 = polyZonotope(center(zono2),generators(zono2),[],[zeros(2,3);eye(3)]);
    
    zono3 = zonotope(rand(2,5)-0.5*ones(2,5));
    pZ3 = polyZonotope(center(zono3),generators(zono3),[],[zeros(5,4);eye(4)]);
    
    % create random third-order tensor
    T{1,1} = rand(2)-0.5*ones(2);
    T{1,2} = rand(2)-0.5*ones(2);
    T{2,1} = rand(2)-0.5*ones(2);
    T{2,2} = rand(2)-0.5*ones(2);
    
    % compute cubic map
    zonoRes = cubMap(zono1,zono2,zono3,T);
    pZres = cubMap(pZ1,pZ2,pZ3,T);
    
    % over approximate the resulting polynomial zonotope with a zonotope
    zonoRes_ = zonotope(pZres);
    
    % test for equality
    Tol = 1e-12;
    
    if any(size(zonoRes_) ~= size(zonoRes_)) || any(abs(center(zonoRes_) - center(zonoRes)) > Tol)
        error('testLongDuration_polyZonotope_cubMap: comparison with zonotope implementation failed!'); 
    else
        G = generators(zonoRes);
        G_ = generators(zonoRes_);
        for i = 1:size(G,2)
            if ~ismembertol(G(:,i)',G_',Tol,'ByRows',true)
                error('testLongDuration_polyZonotope_cubMap: comparison with zonotope implementation failed!'); 
            end
        end
    end
end


% TEST 2-dimensional (mixed multiplication)

for i = 1:3

    % create random polynomial zonotopes
    pZ = cell(3,1);
    
    for j = 1:3
        c = rand(2,1)-0.5*ones(2,1);
        G = rand(2,j+1+5)-0.5*ones(2,j+1+5);
        Grest = rand(2,j+1)-0.5*ones(2,j+1);
        expMat = [eye(j+1), round(rand(j+1,5)*10)];
        pZ{j} = polyZonotope(c,G,Grest,expMat);
    end

    % create random third-order tensor
    T{1,1} = rand(2)-0.5*ones(2);
    T{1,2} = rand(2)-0.5*ones(2);
    T{2,1} = rand(2)-0.5*ones(2);
    T{2,2} = rand(2)-0.5*ones(2);

    % obtain result
    pZres = cubMap(pZ{1},pZ{2},pZ{3},T);

    % determine random points inside the polynomial zonotopes
    N = 20;
    
    points1 = [randPoint(pZ{1},N),randPoint(pZ{1},'all','extreme')];
    points2 = [randPoint(pZ{2},N),randPoint(pZ{2},'all','extreme')];
    points3 = [randPoint(pZ{3},N),randPoint(pZ{3},'all','extreme')];
    
    % calculate the cubic map for all possible point combinations
    pointsRes = zeros(2,size(points1,2)*size(points2,2)*size(points3,2));
    counter = 1;
    
    for j = 1:size(points1,2)
        for k = 1:size(points2,2)
            for h = 1:size(points3,2)
                pointsRes(:,counter) = cubMulPoint(points1(:,j),points2(:,k),points2(:,k),T);
                counter = counter + 1;
            end
        end
    end

%     % visulaize result
%     figure
%     plot(pZres,[1,2],'FaceColor',[.5 .5 .5],'Filled',true,'Order',5,'Splits',30)
%     hold on
%     plot(zonotope(pZres),[1,2],'r');
%     plot(pointsRes(1,:),pointsRes(2,:),'.k');

    % check if the polynomial zonotope obtained by computing the quadratic map
    % encloses all randomly drawn and mapped points
    suc = containsPointSet(pZres,pointsRes,[],30);
   
    if ~suc
        error('testLongDuration_polyZonotope_cubMap: random test mixed multiplication failed!');
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