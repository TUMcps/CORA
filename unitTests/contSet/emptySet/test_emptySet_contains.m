function res = test_emptySet_contains
% test_emptySet_contains - unit test function of contains
%
% Syntax:
%    res = test_emptySet_contains
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

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       05-April-2023
% Last update:   20-January-2025 (AK, added general containment checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init empty set
n = 2;
O = emptySet(n);

% empty
p = double.empty(n,0);
assert(contains(O,p));

% point
p = [2;1];
assert(~contains(O,p));

% zonotope
Z = zonotope(zeros(n,1),eye(n));
assert(~contains(O,Z));

S = emptySet(2);
Sdeg = S;
Sempty = S;

implementedSets = {'capsule','conPolyZono','conZonotope','interval','polytope',...
                    'zonoBundle','zonotope','ellipsoid','taylm',...
                    'polyZonotope','conPolyZono','spectraShadow','levelSet'};

setsNonExact = {};

additionalAlgorithms = {};

additionalAlgorithms_specificSets = {};

checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
