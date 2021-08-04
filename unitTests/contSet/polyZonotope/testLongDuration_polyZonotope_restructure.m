function res = testLongDuration_polyZonotope_restructure
% testLongDuration_polyZonotope_restructure - unit test function for
%    over-approximative polynomial zonotope restructuring
%
% Syntax:  
%    res = testLongDuration_polyZonotope_restructure
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
% Written:      29-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;


%% RANDOM TESTS

% TEST 2-dimensional

methods = {'reduceGirard','reduceMethC','zonotopeGirard'};

for j = 1:length(methods)
    for i = 1:3

        % create random zonotope
        c = rand(2,1)-0.5*ones(2,1);
        G = rand(2,9)-0.5*ones(2,9);
        ind = datasample(1:7,4,'Replace',false);
        G(:,ind) = G(:,ind)./10;
        Grest = (rand(2,3)-0.5*ones(2,3))*2;
        expMat = [eye(4), zeros(4,5)];
        for k = 5:9
           ind1 = 1 + round(rand()*3); 
           ind2 = 1 + round(rand()*3);
           expMat(ind1,k) = round(rand()*5);
           expMat(ind2,k) = round(rand()*5);
        end
        pZ = polyZonotope(c,G,Grest,expMat);

        % restructure the polynomial zonotope
        if strcmp(methods{j},'reduceGirard') && i >= 2
            pZres = restructure(pZ,methods{j},2,3);
        else
            pZres = restructure(pZ,methods{j},2);
        end

        % determine random point and extreme points inside the original polynomial
        % zonotope
        N = 10000;
        points = randPoint(pZ,N);
        pointsExt = randPoint(pZ,'all','extreme');

        points = [pointsExt,points];

        % check if the all points from the original polynomial zonotope are
        % enclosed by the reduced polynomial zonotope
        suc = containsPointSet(pZres,points,[],30);

%         % visualization
%         plot(pZres,[1,2],'FaceColor',[.5 .5 .5],'Filled',true,'Splits',5,'Order',30)
%         hold on
%         plot(zonotope(pZres),[1,2],'r');
%         plot(points(1,:),points(2,:),'.k');

        if ~suc
           error('testLongDuration_polyZonotope_restructure: random test 2D failed!'); 
        end
    end
end


% TEST 4-dimensional

for j = 1:length(methods)
    for i = 1:3

        % create random zonotope
        c = rand(4,1)-0.5*ones(4,1);
        G = rand(4,8)-0.5*ones(4,8);
        ind = datasample(1:6,4,'Replace',false);
        G(:,ind) = G(:,ind)./10;
        Grest = (rand(4,5)-0.5*ones(4,5))*2;
        expMat = [eye(6), round(rand(6,2)*5)];
        pZ = polyZonotope(c,G,Grest,expMat);

        % restructure the polynomial zonotope
        pZres = restructure(pZ,methods{j},1);

        % determine random point and extreme points inside the original polynomial
        % zonotope
        N = 10000;
        points = randPoint(pZ,N);
        pointsExt = randPoint(pZ,'all','extreme');

        points = [pointsExt,points];

        % check if the all points from the original polynomial zonotope are
        % enclosed by the reduced polynomial zonotope
        suc = containsPointSet(pZres,points,[],30);

        if ~suc
           error('testLongDuration_polyZonotope_restructure: random test 4D failed!'); 
        end
    end
end

res = true;

%------------- END OF CODE --------------