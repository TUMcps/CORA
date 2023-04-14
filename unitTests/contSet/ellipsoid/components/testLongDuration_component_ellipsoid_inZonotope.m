function res = testLongDuration_component_ellipsoid_inZonotope
% testLongDuration_component_ellipsoid_inZonotope - unit test function of
%    testLongDuration_ellipsoid_inZonotope
%
% Syntax:  
%    res = testLongDuration_component_ellipsoid_inZonotope
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
% Written:      17-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
for i=2:3
    for j=1:nRuns
        % do not test vertices in high dims as comp. complexity prevents it
        % for high dimensions
        Z = zonotope.generateRandom('Dimension',2);
        E_Zo = ellipsoid(Z,'outer:exact');
        if ~contains(E_Zo,Z)
            res = false;
            break;
        end
        E_Zi = ellipsoid(Z,'inner:exact');
        if ~contains(E_Zo,E_Zi)
            res = false;
            break;
        end
    end
    if ~res
        path = pathFailedTests(mfilename());
        save(path,'Z');
        break;
    end
end
%------------- END OF CODE --------------