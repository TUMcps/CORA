function res = testLongDuration_conPolyZono_reduce
% test_conPolyZono_reduce - unit test function for generator reduction of
%                           constrained polynomial zonotopes
%
% Syntax:  
%    res = test_conPolyZono_reduce()
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
% See also: conPolyZono/reduce

% Author:       Niklas Kochdumper
% Written:      26-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;
    splits = 4;
    
    
    % Random Tests --------------------------------------------------------
    
    % different reduction methods
    methods = {'girard','pca','scott'};
    
    % loop over all test cases
    for i = 1:5
        
        % generate random constrained polynomial zonotope
        cPZ1 = conPolyZono.generateRandom(2,randi([10,20]));
        
        % draw desired reduced order at random
        temp = size(cPZ1.G,2) + size(cPZ1.Grest,2) + size(cPZ1.A,2);
        order = rand() * temp/2;
        
        % compute random points inside the original set
        points = randPoint(cPZ1,10,'extreme');

        % loop over all reduction methods
        for j = 1:length(methods)
        
            % reduce constraints
            cPZ = reduce(cPZ1,methods{j},order);
        
            % check if all points are inside polygon enclosures
            pgon = polygon(cPZ,splits);
            
            if ~in(pgon,points)

                % save variables so that failure can be reproduced
                path = pathFailedTests(mfilename());
                save(path,'cPZ1','order','points');

                error('Random test failed!');
            end
        end
    end
end

%------------- END OF CODE --------------
