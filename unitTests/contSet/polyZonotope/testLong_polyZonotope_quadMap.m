function res = testLong_polyZonotope_quadMap
% testLong_polyZonotope_quadMap - unit test function of quadMap
%
% Syntax:
%    res = testLong_polyZonotope_quadMap
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
% Written:       23-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% TEST 2-dimensional (quadratic multiplication)

for i = 1:3

    % create random zonotope
    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,8)-0.5*ones(2,8);
    GI = rand(2,3)-0.5*ones(2,3);
    E = [eye(3), round(rand(3,5)*10)];
    pZ = polyZonotope(c,G,GI,E);

    % create random matrices
    Q{1} = rand(2)-0.5*ones(2);
    Q{2} = rand(2)-0.5*ones(2);

    % obtain result
    pZres = quadMap(pZ,Q);

    % draw random points from zonotope and approximately compute the
    % quadratic mapping
    N = 10000;
    points = zeros(2,N);

    initPoints = [randPoint(pZ,N),randPoint(pZ,'all','extreme')];

    for j = 1:size(initPoints,2)
        p = initPoints(:,j);
        points(:,j) = quadMapPoint(p,p,Q);
    end

%     % visulaize result
%     plot(pZ,[1,2],'FaceColor',[.5 .5 .5],'Splits',5,'Order',30);
%     hold on
%     plot(initPoints(1,:),initPoints(2,:),'.k');
%     plot(zonotope(pZ),[1,2],'r');
% 
%     figure
%     plot(pZres,[1,2],'FaceColor',[.5 .5 .5],'Splits',5,'Order',30)
%     hold on
%     plot(zonotope(pZres),[1,2],'r');
%     plot(points(1,:),points(2,:),'.k');

    % check if the polynomial zonotope obtained by computing the quadratic map
    % encloses all randomly drawn and mapped points
    if ~containsPointSet(pZres,points,[],30)
        res = false; return
    end  
end


% TEST 2-dimensional (mixed multiplication)

for i = 1:3

    % create tow random zonotopes
    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,8)-0.5*ones(2,8);
    GI = rand(2,3)-0.5*ones(2,3);
    E = [eye(3), round(rand(3,5)*10)];
    pZ1 = polyZonotope(c,G,GI,E);

    c = rand(2,1)-0.5*ones(2,1);
    G = rand(2,6)-0.5*ones(2,6);
    GI = rand(2,2)-0.5*ones(2,2);
    E = [eye(3), round(rand(3,3)*10)];
    pZ2 = polyZonotope(c,G,GI,E);

    % create random matrices
    Q{1} = rand(2)-0.5*ones(2);
    Q{2} = rand(2)-0.5*ones(2);

    % obtain result
    pZres = quadMap(pZ1,pZ2,Q);

    % draw random points from zonotope and approximately compute the
    % quadratic mapping
    N = 100;

    initPoints1 = [randPoint(pZ1,N),randPoint(pZ1,'all','extreme')];
    initPoints2 = [randPoint(pZ2,N),randPoint(pZ2,'all','extreme')];
    
    points = zeros(2,size(initPoints1,2)*size(initPoints2,2));

    counter = 1;
    
    for j = 1:size(initPoints1,2)
        for k = 1:size(initPoints2,2)
            p1 = initPoints1(:,j);
            p2 = initPoints2(:,k);
            points(:,counter) = quadMapPoint(p1,p2,Q);
            counter = counter + 1;
        end
    end

%     % visualize result
%     figure
%     plot(pZres,[1,2],'FaceColor',[.5 .5 .5],'Splits',5,'Order',30)
%     hold on
%     plot(zonotope(pZres),[1,2],'r');
%     plot(points(1,:),points(2,:),'.k');

    % check if the polynomial zonotope obtained by computing the quadratic map
    % encloses all randomly drawn and mapped points
    if ~containsPointSet(pZres,points,[],30)
        res = false; return
    end  
end


% TEST 4-dimensional (quadratic multiplication)

for i = 1:3

    % create random zonotope
    c = rand(4,1)-0.5*ones(4,1);
    G = rand(4,8)-0.5*ones(4,8);
    GI = rand(4,3)-0.5*ones(4,3);
    E = [eye(5), round(rand(5,3)*10)];
    pZ = polyZonotope(c,G,GI,E);

    % create random matrices
    Q{1} = rand(4)-0.5*ones(4);
    Q{2} = rand(4)-0.5*ones(4);
    Q{3} = rand(4)-0.5*ones(4);
    Q{4} = rand(4)-0.5*ones(4);

    % obtain result
    pZres = quadMap(pZ,Q);

    % draw random points from zonotope and approximately compute the
    % quadratic mapping
    N = 10000;
    points = zeros(4,N);

    initPoints = [randPoint(pZ,N),randPoint(pZ,'all','extreme')];

    for j = 1:size(initPoints,2)
        p = initPoints(:,j);
        points(:,j) = quadMapPoint(p,p,Q);
    end

    % check if the polynomial zonotope obtained by computing the quadratic map
    % encloses all randomly drawn and mapped points
    if ~containsPointSet(pZres,points,[],5)
        res = false; return
    end  
end

% ------------------------------ END OF CODE ------------------------------
