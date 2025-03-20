function res = testLong_conPolyZono_conZonotope
% testLong_conPolyZono_conZonotope - unit test function for 
%    constrained zonotope enclosure of constrained polynomial zonotopes
%
% Syntax:
%    testLong_conPolyZono_conZonotope
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
% See also: conPolyZono/conZonotope

% Authors:       Niklas Kochdumper
% Written:       26-January-2021
% Last update:   05-March-2025 (TL, reduced runtime)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
tol = 1e-5;

% define range bounding methods that are tested
methods = {'extend','linearize','all'};

% loop over all test cases
for i = 1:5
    
    % generate random constrained polynomial zonotope
    cPZ = conPolyZono.generateRandom();
    
    % get random points
    points = randPoint(cPZ,10,'extreme');
    
    % loop over all methods
    for j = 1:length(methods)
        
        % compute zonotope enclosure
        cZ = conZonotope(cPZ,methods{j});

        % reduce before containment check
        cZ = reduce(cZ,'girard',2);
        
        % check for correctness
        assertLoop(contains(cZ,points,'exact',tol),i,j)
    end
end

% ------------------------------ END OF CODE ------------------------------
