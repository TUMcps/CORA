function res = testFlaky_ellipsoid_isIntersecting
% testFlaky_ellipsoid_isIntersecting - unit test function of
%    isIntersecting
%
% Syntax:
%    res = testFlaky_ellipsoid_isIntersecting
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
% See also: none

% Authors:       Victor Gassmann
% Written:       18-March-2021
% Last update:   ---
% Last revision: 22-September-2024 (MW, integrate components)

% ------------------------------ BEGIN CODE -------------------------------

nRuns = 2;

% degeneracy
bools = [false,true];

% smaller dims since halfspaces and vertices are involved
for i=5:5:10
    for j=1:nRuns
        for k=1:2 
            % init random ellipsoid, zonotope, and polyZonotope
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            Z = zonotope.generateRandom('Dimension',i);
            pZ = polyZonotope.generateRandom('Dimension',i);

            % sample random point from ellipsoid
            s = randPoint(E,1);

            % init zonotope and polyZonotope containing that random point
            % so that sets always intersect
            Z1 = zonotope(s,Z.G);
            pZ1 = polyZonotope(s,pZ.G, pZ.GI, pZ.E, pZ.id);
            
            % check for intersection
            assertLoop(isIntersecting(E,Z1,'approx'),i,j,k)
            assertLoop(isIntersecting(E,pZ1,'approx'),i,j,k)
        end
    end
end

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
