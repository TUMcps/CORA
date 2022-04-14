function res = testLongDuration_conPolyZono_and
% test_conPolyZono_and - unit test function for intersection of constrained
%                        polynomial zonotopes
%
% Syntax:  
%    res = test_conPolyZono_and()
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
% See also: conPolyZono/and

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
        cPZ1 = conPolyZono.generateRandom(2,[],[],[],0);
        
        % loop over all other set representations
        for j = 1:length(sets)
            
            found = 0;
            
            % loop unit two sets which really intersect are obtained
            for h = 1:3
            
                % generate random object of the current set representation
                if strcmp(sets{j},'conPolyZono')
                    temp = conPolyZono.generateRandom(2,[],[],[],0);
                elseif strcmp(sets{j},'polyZonotope')
                    temp = polyZonotope.generateRandom(2,[],[],0);
                else
                    str = ['temp = ',sets{j},'.generateRandom(2);']; 
                    eval(str);
                end
                
                cPZ2 = conPolyZono(temp);

                % compute intersection
                cPZ = cPZ1 & temp;

                % get random points
                try
                   points = randPoint(cPZ,10,'extreme');
                catch
                   continue 
                end
                
                found = 1;
                break;
            end
            
            % check if two intersecting sets could be generated
            if found == 0
               continue; 
            end
            
            % compute polytope enclosure of the two conPolyZono objects
            pgon1 = polygon(cPZ1,splits);
            pgon2 = polygon(cPZ2,splits);
            
            pgon = pgon1 & pgon2;
            
            % check if all points are inside the intersection of the
            % polygon enclosures
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