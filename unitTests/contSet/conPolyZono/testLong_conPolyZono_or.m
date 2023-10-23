function res = testLong_conPolyZono_or
% testLong_conPolyZono_or - unit test function for 
%    union of constrained polynomial zonotopes
%
% Syntax:
%    res = testLong_conPolyZono_or()
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
% See also: conPolyZono/or

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
        str = ['temp = ',sets{i},'.generateRandom(''Dimension'',2);']; 
        eval(str);
        cPZ2 = conPolyZono(temp);

        % compute union
        cPZ = cPZ1 | temp;

        % get random points inside the two conPolyZono objects
        points = [randPoint(cPZ1,5,'extreme'), ...
                  randPoint(cPZ2,5,'extreme')];
        
        % check if all points are inside polygon enclosures
        pgon = polygon(cPZ,splits);
        
        if ~contains(pgon,points)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
