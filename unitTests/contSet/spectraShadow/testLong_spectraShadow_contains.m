function res = testLong_spectraShadow_contains
% testLong_spectraShadow_contains - unit test function of contains
%
% Syntax:
%    res = testLong_spectraShadow_contains
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
% Written:       14-August-2023
% Last update:   21-January-2025 (AK, added general containment checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

runs = 5;

% Simple test: Spectrahedral shadow should contain its center
for i=1:runs
    SpS = spectraShadow.generateRandom();
    c = center(SpS);
    assert(contains(SpS,c,'exact',1e-4))
end

% The test of whether randPoint-points are contained in S will be done in
% testLong_spectraShadow_randPoint.m, so no need to do it twice

% Check all containment combinations
I = interval([-1;-1], [1;1]);
Ideg = interval([-1;0],[1;0]);

S = spectraShadow(I);
Sdeg = spectraShadow(Ideg);
Sempty = spectraShadow.empty(2);

implementedSets = {'capsule','conZonotope','interval','polytope',...
                    'zonoBundle','zonotope','ellipsoid','taylm',...
                    'polyZonotope','conPolyZono','spectraShadow'};

setsNonExact = {'capsule','ellipsoid','taylm','polyZonotope','conPolyZono',...
                'spectraShadow'};

additionalAlgorithms = {'sampling'};
additionalAlgorithms_specificSets = {'conZonotope','interval','polytope',...
                                       'zonoBundle', 'zonotope'}; % for 'sampling'

checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

res = true;

% ------------------------------ END OF CODE ------------------------------
