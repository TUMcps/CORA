function res = testLong_interval_contains
% testLong_interval_contains - unit test function of contains
%
% Syntax:
%    res = testLong_interval_contains
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

% Authors:       Adrian Kulmburg
% Written:       20-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

S = 2*interval([-1;-1], [1;1]);
Sdeg = 2*interval([-1;0], [1;0]);
Sempty = interval.empty(2);

implementedSets = {'capsule','conPolyZono','conZonotope','interval','polytope',...
                    'zonoBundle','zonotope','ellipsoid','taylm',...
                    'polyZonotope','conPolyZono','spectraShadow'};

setsNonExact = {'taylm','conPolyZono','polyZonotope'};

additionalAlgorithms = {};

additionalAlgorithms_specificSets = {};

checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
