function res = testLong_polyZonotope_zonotope
% testLong_polyZonotope_zonotope - unit test function for zonotope
%    over-approximation of a polynomial zonotope
%
% Syntax:
%    res = testLong_polyZonotope_zonotope
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

    % calculate zonotope over-approximation
    Z = zonotope(pZ);
    Z = halfspace(Z);

    % determine random point and extreme points inside the original polynomial
    % zonotope
    N = 10000;
    points = randPoint(pZ,N);
    pointsExt = randPoint(pZ,'all','extreme');

    points = [pointsExt,points];

    % check if the all points from the original polynomial zonotope are
    % enclosed by the reduced polynomial zonotope
    if ~contains(Z, points)
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

    % calculate zonotope over-approximation
    Z = zonotope(pZ);
    Z = halfspace(Z);

    % determine random point and extreme points inside the original
    % polynomial zonotope
    N = 10000;
    points = randPoint(pZ,N);
    pointsExt = randPoint(pZ,'all','extreme');

    points = [pointsExt,points];

    % check if the all points from the original polynomial zonotope are
    % enclosed by the reduced polynomial zonotope
    if ~contains(Z, points)
        res = false; return
    end
    
end

% ------------------------------ END OF CODE ------------------------------
