function res = test_fullspace_contains
% test_fullspace_contains - unit test function of contains
%
% Syntax:
%    res = test_fullspace_contains
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

% init fullspace
n = 2;
fs = fullspace(n);

% init point
p = [2;1];
assert(contains(fs,p));

% init zonotope
Z = zonotope([2;1],eye(n));
assert(contains(fs,Z));

% init ellipsoid
E = ellipsoid(eye(n),ones(n,1));
assert(contains(fs,E));

% empty set
O = emptySet(n);
assert(contains(fs,O));

S = fullspace(2);
% Fullspace has no degenerate/empty representation, thus both need to be
% replaced
Sdeg = emptySet(2);
Sempty = emptySet(2);

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
