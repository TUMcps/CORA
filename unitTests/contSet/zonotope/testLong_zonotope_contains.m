function res = testLong_zonotope_contains
% testLong_zonotope_contains - unit test function of contains
%
% Syntax:
%    res = testLong_zonotope_contains
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

% Authors:       Matthias Althoff, Adrian Kulmburg
% Written:       26-July-2016
% Last update:   14-September-2019
%                01-July-2021 (AK, more tests, and merged with testLong_zonotope_containsPoint)
%                21-January-2025 (AK, added general containment checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% zonotope x parallelotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);
P1 = zonotope([-3.8, -4, 3; 1.2, 3, -4]);
P2 = zonotope([-3.8, -8, 2; 1.2, 10, -10]);

% containment test with all methods ---
assert(~contains(P1,Z1));
assert(contains(P2,Z1));

% exact
assert(~contains(P1,Z1,'exact'));
assert(contains(P2,Z1,'exact'));

assert(~contains(P1,Z1,'exact:venum'));
assert(contains(P2,Z1,'exact:venum'));

assert(~contains(P1,Z1,'exact:polymax'));
assert(contains(P2,Z1,'exact:polymax'));

% opt
assert(~contains(P1,Z1,'opt',0,200));
assert(contains(P2,Z1,'opt',0,200));

% approx
assert(~contains(P1,Z1,'approx'));
assert(contains(P2,Z1,'approx'));

assert(~contains(P1,Z1,'approx:st'));
assert(contains(P2,Z1,'approx:st'));

assert(~contains(P1,Z1,'approx:stDual'));
assert(contains(P2,Z1,'approx:stDual'));

% sampling
assert(~contains(P1,Z1,'sampling',0,200));
assert(contains(P2,Z1,'sampling',0,200));

assert(~contains(P1,Z1,'sampling:primal',0,200));
assert(contains(P2,Z1,'sampling:primal',0,200));

assert(~contains(P1,Z1,'sampling:dual',0,200));
assert(contains(P2,Z1,'sampling:dual',0,200));


% zonotope x point
n = 2; nrGen = 5;
Z = zonotope(zeros(n,1),rand(n,nrGen));

% point inside: center of zonotope
p_inside = center(Z);
% point outside: add all generators (with max of rand -> 1)
p_outside = nrGen*ones(n,1);
% array of points
num = 10;
p_array = nrGen*(ones(n,num)+rand(n,num));

% check if correct results for containment
assert(contains(Z, p_inside));
assert(~contains(Z, p_outside));
assert(~all(contains(Z,p_array)));

% Check all containment combinations
I = interval([-1;-1],[1;1]);
Ideg = interval([-1;0], [1;0]);

S = zonotope(I);
Sdeg = zonotope(Ideg);
Sempty = zonotope.empty(2);

implementedSets = {'capsule','conPolyZono','conZonotope','interval','polytope',...
                    'zonoBundle','zonotope','ellipsoid','taylm',...
                    'polyZonotope','conPolyZono','spectraShadow'};

setsNonExact = {'taylm','conPolyZono','polyZonotope'}; 
% 'spectraShadow' is generally not exact, but for the specific zonotope (an interval), it is.
% However, checkAllContainments expects it that way...

additionalAlgorithms = {'exact:venum', 'exact:polymax', 'approx:st', 'approx:stDual', 'sampling', 'sampling:primal', 'sampling:dual'};

additionalAlgorithms_specificSets = {{'conZonotope','interval', ...
                                      'polytope', 'zonoBundle', ...
                                      'zonotope'},... % exact:venum
                                     {'ellipsoid','conZonotope', ...
                                      'interval','polytope',...
                                      'zonoBundle','zonotope'},... % exact:polymax
                                     {'zonotope'},... % approx:st
                                     {'zonotope'},... % approx:stDual
                                     {'zonotope'},... % sampling
                                     {'zonotope'},... % sampling:primal
                                     {'zonotope'},... % sampling:dual
    };

checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);


% all good
res = true;

% ------------------------------ END OF CODE ------------------------------
