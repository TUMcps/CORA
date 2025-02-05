function res = testLong_conPolyZono_contains
% testLong_conPolyZono_contains - unit test function for containment
%    checks of constrained polynomial zonotopes
%
% Syntax:
%    res = testLong_conPolyZono_contains
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

% Authors:       Adrian Kulmburg
% Written:       17-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Check all containment combinations
S = conPolyZono(interval([-1;-1], [1;1]));
Sdeg = conPolyZono(interval([-1;0], [1;0]));
Sempty = conPolyZono.empty(2);

implementedSets = {'interval', 'zonotope', 'polytope', 'zonoBundle',...
                    'conZonotope', 'taylm', 'polyZonotope', 'conPolyZono'};

setsNonExact = [implementedSets {'point'}];

additionalAlgorithms = {};

additionalAlgorithms_specificSets = {};

checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

res = true;

% ------------------------------ END OF CODE ------------------------------
