function res = testLong_conPolyZono_convHull
% testLong_conPolyZono_convHull - unit test function for the 
%    convex hull of constrained polynomial zonotopes
%
% Syntax:
%    res = testLong_conPolyZono_convHull()
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
splits = 4;

% Random Tests --------------------------------------------------------

% define set representations that are tested
sets = {'conPolyZono','polyZonotope','zonotope','conZonotope', ...
        'ellipsoid','capsule'};

% loop over all test cases
for i = 1:3
    
    % generate random constrained polynomial zonotope
    cPZ1 = conPolyZono.generateRandom('Dimension',2);
    
    % loop over all other set representations
    for j = 1:length(sets)
        
        % generate random object of the current set representation
        eval(['temp = ',sets{i},'.generateRandom(''Dimension'',2);']);
        cPZ2 = conPolyZono(temp);

        % compute convex hull
        cPZ = convHull(cPZ1,temp);

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
        
        if ~contains(pgon,points)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
