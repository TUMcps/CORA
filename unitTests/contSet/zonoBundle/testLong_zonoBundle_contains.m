function res = testLong_zonoBundle_contains
% testLong_zonoBundle_contains - unit test function of contains
%
% Syntax:
%    res = testLong_zonoBundle_contains
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
% Written:       21-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Check all containment combinations
I = interval([-1;-1],[1;1]);
Ideg = interval([-1;0], [1;0]);

S = zonoBundle(I);
Sdeg = zonoBundle(Ideg);
Sempty = zonoBundle.empty(2);

setsExact = {'capsule','conZonotope','interval','polytope',...
             'zonoBundle','zonotope','ellipsoid',...
             'spectraShadow'};

setsNonExact = {};

implementedSets = [setsExact setsNonExact];

additionalAlgorithms = {'exact:zonotope', 'exact:polytope'};

additionalAlgorithms_specificSets = {setsExact,... % exact:zonotope
                                     setsExact,... % exact:polytope
    };

checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

res = true;

% ------------------------------ END OF CODE ------------------------------
