function res = test_conPolyZono_convHull
% test_conPolyZono_convHull - unit test function for the 
%    convex hull of constrained polynomial zonotopes
%
% Syntax:
%    res = test_conPolyZono_convHull()
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
% See also: conPolyZono/convHull

% Authors:       Niklas Kochdumper
% Written:       03-February-2021
% Last update:   22-May-2025 (TL, converted to short test)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
splits = 4;

% Init sets ---------------------------------------------------------------

% cPZ
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];
cPZ1 = conPolyZono(c,G,E,A,b,EC);

% Z
c = [1;1];
G = [1 1 1; 1 -1 0];
Z = zonotope(c,G);

% define set representations that are tested
sets = {cPZ1,Z};

% loop over all other set representations
for j = 1:length(sets)
    
    % generate random object of the current set representation
    S = sets{j};
    cPZ2 = conPolyZono(S);

    % compute convex hull
    cPZ = convHull(cPZ1,S);

    % get random points inside the two conPolyZono objects
    N1 = 10;
    points1 = [randPoint(cPZ1,N1/2,'extreme'), ...
               randPoint(cPZ2,N1/2,'extreme')];
    
    % compute random combinations of the points
    N2 = 100;
    points2 = zeros(dim(cPZ),N2);
    
    for k = 1:N2
        ind1 = randi([1,N1]);
        ind2 = randi([1,N1]);
        points2(:,k) = points1(:,ind1) + ...
                       rand()*(points1(:,ind2) - points1(:,ind1));
    end
    
    points = [points1,points2];
    
    % check if all points are inside polygon enclosures
    pgon = polygon(cPZ,splits);
    
    assertLoop(contains(pgon,points),j)
end

end

% ------------------------------ END OF CODE ------------------------------
