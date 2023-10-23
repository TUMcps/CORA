function res = testLong_conPolyZono_linComb
% testLong_conPolyZono_linComb - unit test function for the 
%    linear combination of constrained polynomial zonotopes
%
% Syntax:
%    res = testLong_conPolyZono_linComb()
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
% See also: conPolyZono/linComb

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
for i = 1:2
    
    % generate random constrained polynomial zonotope
    cPZ1 = conPolyZono.generateRandom('Dimension',2);
    
    % loop over all other set representations
    for j = 1:length(sets)
        
        % generate random object of the current set representation
        eval(['temp = ',sets{i},'.generateRandom(''Dimension'',2);']);
        cPZ2 = conPolyZono(temp);

        % compute linear combination
        cPZ = linComb(cPZ1,temp);

        % get random points inside the two conPolyZono objects
        N1 = 10;
        points1 = randPoint(cPZ1,N1,'extreme');
        points2 = randPoint(cPZ2,N1,'extreme');
        
        % compute random combinations of the points
        N2 = 100;
        points3 = zeros(dim(cPZ),N2);
        
        for k = 1:N2
            ind1 = randi([1,N1]);
            ind2 = randi([1,N1]);
            points2(:,k) = points1(:,ind1) + ...
                           rand()*(points2(:,ind2) - points1(:,ind1));
        end
        
        % concatenate points
        points = [points1,points2,points3];
        
        % check if all points are inside polygon enclosures
        pgon = polygon(cPZ,splits);
        
        if ~contains(pgon,points)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
