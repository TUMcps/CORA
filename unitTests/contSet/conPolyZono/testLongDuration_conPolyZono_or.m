function res = testLongDuration_conPolyZono_or
% test_conPolyZono_or - unit test function for union of constrained
%                       polynomial zonotopes
%
% Syntax:  
%    res = test_conPolyZono_or()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conPolyZono/or

% Author:       Niklas Kochdumper
% Written:      03-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;
    splits = 4;
    
    % Random Tests --------------------------------------------------------
    
    % define set representations that are tested
    sets = {'conPolyZono','polyZonotope','zonotope','conZonotope', ...
            'ellipsoid','capsule'};
    
    % loop over all test cases
    for i = 1:2
        
        % generate random constrained polynomial zonotope
        cPZ1 = conPolyZono.generateRandom(2);
        
        % loop over all other set representations
        for j = 1:length(sets)
            
            % generate random object of the current set representation
            str = ['temp = ',sets{i},'.generateRandom(2);']; 
            eval(str);
            cPZ2 = conPolyZono(temp);

            % compute union
            cPZ = cPZ1 | temp;

            % get random points inside the two conPolyZono objects
            points = [randPoint(cPZ1,5,'extreme'), ...
                      randPoint(cPZ2,5,'extreme')];
            
            % check if all points are inside polygon enclosures
            pgon = polygon(cPZ,splits);
            
            if ~in(pgon,points)
                
                % save variables so that failure can be reproduced
                path = pathFailedTests(mfilename());
                save(path,'cPZ1','cPZ2','points');
                
                error('Random test failed!');
            end
        end
    end
end

%------------- END OF CODE --------------
