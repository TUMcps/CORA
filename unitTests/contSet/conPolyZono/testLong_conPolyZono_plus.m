function res = testLong_conPolyZono_plus
% testLong_conPolyZono_plus - unit test function for 
%    Minkowski addition of constrained polynomial zonotopes
%
% Syntax:
%    res = testLong_conPolyZono_plus()
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
% See also: conPolyZono/plus

% Authors:       Niklas Kochdumper
% Written:       03-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
splits = 4;

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

        % compute union
        cPZ = cPZ1 + temp;

        % get random points inside the two conPolyZono objects
        N = 10;
        points1 = randPoint(cPZ1,N,'extreme');
        points2 = randPoint(cPZ2,N,'extreme');
        
        % add the points
        points = zeros(dim(cPZ),N^2);
        cnt = 1;
        
        for k = 1:N
            for l = 1:N
                points(:,cnt) = points1(:,k) + points2(:,l);
                cnt = cnt + 1;
            end
        end
        
        % check if all points are inside polygon enclosures
        pgon = polygon(cPZ,splits);
        
        if ~contains(pgon,points)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
