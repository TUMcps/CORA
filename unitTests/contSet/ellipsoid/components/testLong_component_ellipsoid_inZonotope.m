function res = testLong_component_ellipsoid_inZonotope
% testLong_component_ellipsoid_inZonotope - unit test function of
%    testLong_ellipsoid_inZonotope
%
% Syntax:
%    res = testLong_component_ellipsoid_inZonotope
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
% Written:       17-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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
            return;
        end
        E_Zi = ellipsoid(Z,'inner:exact');
        if ~contains(E_Zo,E_Zi)
            res = false;
            return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
