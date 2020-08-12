function res = test_polyZonotope_quadMap
% test_polyZonotope_quadMap - unit test function of quadMap
%
% Syntax:  
%    res = test_polyZonotope_quadMap
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
% Written:      23-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

%% ANALYTICAL TESTS

% create polynomial matrices
c = [1;2];
G = [1 -2 1;2 3 -1];
Grest = [0;0];
expMat = [1 0 2; 0 1 1];
pZ = polyZonotope(c,G,Grest,expMat);

% create matrices of the quadratic map
Q{1} = [1 2;-1 2];
Q{2} = [-3 0;1 1];

% calculate quadratic map
pZres = quadMap(pZ,Q);

% define ground truth
G = [22 19 -5 11 19 -5 16 -11 2; 6 23 -9 3 23 -9 -9 11 -3];
c = [11;3];
expMat = [1 0 2 2 1 3 0 2 4; 0 1 1 0 1 1 2 2 2];

% check for correctness
if ~all(pZres.c == c)
    error('test_polyZonotope_quadMap: analytical test failed!');
end

for i = 1:size(expMat,2)
    
    ind = ismember(pZres.expMat',expMat(:,i)','rows');  
    ind_ = find(ind > 0);
    
    if isempty(ind_)
        error('test_polyZonotope_quadMap: analytical test failed!');        
    elseif ~all(pZres.G(:,ind_(1)) == G(:,i))
        error('test_polyZonotope_quadMap: analytical test failed!');
    end
    
end






%% RANDOM TESTS

% TEST 2-dimensional (quadratic multiplication)

for i = 1:3

    % create random zonotope
    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,8)-0.5*ones(2,8);
    Grest = rand(2,3)-0.5*ones(2,3);
    expMat = [eye(3), round(rand(3,5)*10)];
    pZ = polyZonotope(c,G,Grest,expMat);

    % create random matrices
    Q{1} = rand(2)-0.5*ones(2);
    Q{2} = rand(2)-0.5*ones(2);

    % obtain result
    pZres = quadMap(pZ,Q);

    % draw random points from zontope and approximately compute the
    % quadratic mapping
    N = 10000;
    points = zeros(2,N);

    initPoints = [pointSet(pZ,N),pointSetExtreme(pZ)];

    for j = 1:size(initPoints,2)
        p = initPoints(:,j);
        value = zeros(2,1);
        value(1) = p' * Q{1} * p;
        value(2) = p' * Q{2} * p;
        points(:,j) = value;
    end

%     % visulaize result
%     plot(pZ,[1,2],'FaceColor',[.5 .5 .5],'Filled',true,'Splits',5,'Order',30);
%     hold on
%     plot(initPoints(1,:),initPoints(2,:),'.k');
%     plot(zonotope(pZ),[1,2],'r');
% 
%     figure
%     plot(pZres,[1,2],'FaceColor',[.5 .5 .5],'Filled',true,'Splits',5,'Order',30)
%     hold on
%     plot(zonotope(pZres),[1,2],'r');
%     plot(points(1,:),points(2,:),'.k');

    % check if the polynomial zonotope obtained by computing the quadratic map
    % encloses all randomly drawn and mapped points
    suc = containsPointSet(pZres,points,[],30);
   
    if ~suc
        error('test_polyZonotope_quadMap: random test 2D failed!');
    end  
end



% TEST 2-dimensional (mixed multiplication)

for i = 1:3

    % create tow random zonotopes
    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,8)-0.5*ones(2,8);
    Grest = rand(2,3)-0.5*ones(2,3);
    expMat = [eye(3), round(rand(3,5)*10)];
    pZ1 = polyZonotope(c,G,Grest,expMat);

    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,6)-0.5*ones(2,6);
    Grest = rand(2,2)-0.5*ones(2,2);
    expMat = [eye(3), round(rand(3,3)*10)];
    pZ2 = polyZonotope(c,G,Grest,expMat);

    % create random matrices
    Q{1} = rand(2)-0.5*ones(2);
    Q{2} = rand(2)-0.5*ones(2);

    % obtain result
    pZres = quadMap(pZ1,pZ2,Q);

    % draw random points from zontope and approximately compute the
    % quadratic mapping
    N = 100;

    initPoints1 = [pointSet(pZ1,N),pointSetExtreme(pZ1)];
    initPoints2 = [pointSet(pZ2,N),pointSetExtreme(pZ2)];
    
    points = zeros(2,size(initPoints1,2)*size(initPoints2,2));

    counter = 1;
    
    for j = 1:size(initPoints1,2)
        for k = 1:size(initPoints2,2)
            p1 = initPoints1(:,j);
            p2 = initPoints2(:,k);
            value = zeros(2,1);
            value(1) = p1' * Q{1} * p2;
            value(2) = p1' * Q{2} * p2;
            points(:,counter) = value;
            counter = counter + 1;
        end
    end

%     % visualize result
%     figure
%     plot(pZres,[1,2],'FaceColor',[.5 .5 .5],'Filled',true,'Splits',5,'Order',30)
%     hold on
%     plot(zonotope(pZres),[1,2],'r');
%     plot(points(1,:),points(2,:),'.k');

    % check if the polynomial zonotope obtained by computing the quadratic map
    % encloses all randomly drawn and mapped points
    suc = containsPointSet(pZres,points,[],30);
   
    if ~suc
        error('test_polyZonotope_quadMap: random test 2D (mixed multiplication) failed!');
    end  
end



% TEST 4-dimensional (quadratic multiplication)

for i = 1:3

    % create random zonotope
    c = rand(4,1)-0.5*ones(4,1);
    G = rand(4,8)-0.5*ones(4,8);
    Grest = rand(4,3)-0.5*ones(4,3);
    expMat = [eye(5), round(rand(5,3)*10)];
    pZ = polyZonotope(c,G,Grest,expMat);

    % create random matrices
    Q{1} = rand(4)-0.5*ones(4);
    Q{2} = rand(4)-0.5*ones(4);
    Q{3} = rand(4)-0.5*ones(4);
    Q{4} = rand(4)-0.5*ones(4);

    % obtain result
    pZres = quadMap(pZ,Q);

    % draw random points from zontope and approximately compute the
    % quadratic mapping
    N = 10000;
    points = zeros(4,N);

    initPoints = [pointSet(pZ,N),pointSetExtreme(pZ)];

    for j = 1:size(initPoints,2)
        p = initPoints(:,j);
        value = zeros(4,1);
        value(1) = p' * Q{1} * p;
        value(2) = p' * Q{2} * p;
        value(3) = p' * Q{3} * p;
        value(4) = p' * Q{4} * p;
        points(:,j) = value;
    end

    % check if the polynomial zonotope obtained by computing the quadratic map
    % encloses all randomly drawn and mapped points
    suc = containsPointSet(pZres,points,[],5);
   
    if ~suc
        error('test_polyZonotope_quadMap: random test 4D failed!');
    end  
end
    

res = true;

%------------- END OF CODE --------------