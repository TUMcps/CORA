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

% Author:       Victor Gassmann
% Written:      18-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
bools = [false,true];
% smaller dims since halfspaces and vertices are involved
for i=5:5:10
    for j=1:nRuns
        for k=1:2 
            %%% generate all variables necessary to replicate results
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            Z = zonotope.generateRandom('Dimension',i);
            pZ = polyZonotope.generateRandom('Dimension',i);
            s = randPoint(E,1);
            %%%
            Z1 = zonotope(s,generators(Z));
            pZ1 = polyZonotope(s,pZ.G,pZ.Grest,pZ.expMat,pZ.id);
            
            % intersect for sure
            if ~isIntersecting(E,Z1,'approx') || ~isIntersecting(E,pZ1,'approx')
                res = false;
                return;
            end
        end
    end
end
%------------- END OF CODE --------------