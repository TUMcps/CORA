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
% Last update:   ---
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
    points = randPoint(cPZ,50,'extreme');
    
    % loop over all methods
    for j = 1:length(methods)
        
        % compute zonotope enclosure
        cZ = conZonotope(cPZ,methods{j});
        
        % check for correctness
        if ~contains(cZ,points,'exact',tol)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
