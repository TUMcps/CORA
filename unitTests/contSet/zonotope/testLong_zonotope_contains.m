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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% zonotope x parallelotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);
P1 = zonotope([-3.8, -4, 3; 1.2, 3, -4]);
P2 = zonotope([-3.8, -8, 2; 1.2, 10, -10]);

% containment test with all methods
assert(~contains(P1,Z1));
assert(contains(P2,Z1));

assert(~contains(P1,Z1,'venum'));
assert(contains(P2,Z1,'venum'));

assert(~contains(P1,Z1,'polymax'));
assert(contains(P2,Z1,'polymax'));

assert(~contains(P1,Z1,'opt',0,200));
assert(contains(P2,Z1,'opt',0,200));

assert(~contains(P1,Z1,'st'));
assert(contains(P2,Z1,'st'));


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


% all good
res = true;

% ------------------------------ END OF CODE ------------------------------
