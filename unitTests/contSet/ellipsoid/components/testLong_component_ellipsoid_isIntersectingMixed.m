function res = testLong_component_ellipsoid_isIntersectingMixed
% testLong_component_ellipsoid_isIntersectingMixed - unit test 
%    function of ellipsoid/isIntersectingMixed
%
% Syntax:
%    res = testLong_component_ellipsoid_isIntersectingMixed
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
% See also: -

% Authors:       Victor Gassmann
% Written:       18-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
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
            if ~isIntersecting(E,Z1,'approx') || ~isIntersecting(E,pZ1,'approx')
                res = false;
                return;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
