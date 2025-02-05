function res = testLong_conZonotope_contains
% testLong_conZonotope_contains - unit test function for containment check
%
% Syntax:
%    res = testLong_conZonotope_contains
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
I = interval([-1;-1],[1;1]);
Ideg = interval([-1;0], [1;0]);

S = conZonotope(I);
Sdeg = conZonotope(Ideg);
Sempty = conZonotope.empty(2);

implementedSets = {'capsule','conPolyZono','conZonotope','interval','polytope',...
                    'zonoBundle','zonotope','ellipsoid','taylm',...
                    'polyZonotope','conPolyZono','spectraShadow'};

setsNonExact = {'taylm','conPolyZono','polyZonotope',...
                'spectraShadow'};

additionalAlgorithms = {'exact:venum', 'exact:polymax', 'approx:st', 'sampling', 'sampling:primal', 'sampling:dual'};

additionalAlgorithms_specificSets = {{'conZonotope','interval', ...
                                      'polytope', 'zonoBundle', ...
                                      'zonotope'},... % exact:venum
                                     {'ellipsoid','conZonotope', ...
                                      'interval','polytope',...
                                      'zonoBundle','zonotope'},... % exact:polymax
                                     {'interval', 'polytope', ...
                                      'zonoBundle', 'conZonotope', ...
                                      'zonotope'},... % approx:st
                                     {'conZonotope'},... % sampling
                                     {'conZonotope'},... % sampling:primal
                                     {'conZonotope'},... % sampling:dual
    };

checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

res = true;


% ------------------------------ END OF CODE ------------------------------
