function res = testLongDuration_conPolyZono_conZonotope
% test_conPolyZono_conZonotope - unit test function for constrained 
%                                zonotope enclosure of constrained 
%                                polynomial zonotopes
%
% Syntax:  
%    test_conPolyZono_conZonotope
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
% See also: conPolyZono/conZonotope

% Author:       Niklas Kochdumper
% Written:      26-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;
    tol = 1e-5;
    
    % Random Tests --------------------------------------------------------
    
    % define range bounding methods that are tested
    methods = {'extend','linearize','all'};
    
    % loop over all test cases
    for i = 1:5
        
        % generate random constrained polynomial zonotope
        cPZ = conPolyZono.generateRandom();
        
        % get random points
        points = randPoint(cPZ,50,'extreme');
        
        % loop over all methods
        for j = 1:length(methods)
            
            % compute zonotope enclosure
            cZ = conZonotope(cPZ,methods{j});
            
            % check for correctness
            if ~in(cZ,points,'exact',tol)
                
                % save variables so that failure can be reproduced
                path = pathFailedTests(mfilename());
                save(path,'cPZ','points');
                
                error('Random test failed!');
            end
        end
    end
end

%------------- END OF CODE --------------
