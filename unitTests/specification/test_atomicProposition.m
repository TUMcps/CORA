function res = test_atomicProposition
% test_atomicProposition - unit test of atomicProposition class
%
% Syntax:
%    res = test_atomicProposition
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Benedikt Seidl
% Written:       09-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% atomic propositions
ap1 = atomicProposition(polytope([0 1 0], 3));
ap2 = atomicProposition(polytope([-1 0 0], -2));
ap3 = atomicProposition(polytope([-1 0 -1], -9));
ap4 = atomicProposition(ellipsoid([13 7; 7 5], [1; 2]), [1 3]);

% sets
set1 = zonotope([1; 1; 1], [1 1 1; 1 -1 0; 0 1 1]);
set2 = zonotope([-7; 1; 4], [1 1 1; 1 -1 0; 0 1 1]);
set3 = zonotope([1; 1; 2], [0.5 1 0.2; 0.2 -1 0; 0 1 0.2]);

% test
assert(ap1.canBeTrue(set1,1));
assert(~ap1.canBeFalse(set1,1));

assert(ap2.canBeTrue(set1,1));
assert(ap2.canBeFalse(set1,1));

assert(~ap3.canBeTrue(set1,1));
assert(ap3.canBeFalse(set1,1));

assert(ap4.canBeTrue(set1,1));
assert(ap4.canBeFalse(set1,1));

assert(~ap4.canBeTrue(set2,1));
assert(ap4.canBeFalse(set2,1));

assert(ap4.canBeTrue(set3,1));
assert(~ap4.canBeFalse(set3,1));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
