function res = test_conPolyZono_conZonotope
% test_conPolyZono_conZonotope - unit test function for 
%    constrained zonotope enclosure of constrained polynomial zonotopes
%
% Syntax:
%    test_conPolyZono_conZonotope
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
%                22-May-2025 (TL, converted to short test)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-5;

% define range bounding methods that are tested
methods = {'extend','linearize'};
    
% init constrained polynomial zonotope
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];
cPZ = conPolyZono(c,G,E,A,b,EC);

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

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
