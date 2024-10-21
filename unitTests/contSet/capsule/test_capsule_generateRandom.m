function res = test_capsule_generateRandom
% test_capsule_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_capsule_generateRandom
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

% Authors:       Mark Wetzlinger
% Written:       27-September-2019
% Last update:   19-May-2022 (name-value pair syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty call
C = capsule.generateRandom();

% values for tests
n = 3;
c = [2;1;-1];
r = 3.5;

% only dimension
C = capsule.generateRandom('Dimension',n);
assert(dim(C) == n);

% only center
C = capsule.generateRandom('Center',c);
assert(compareMatrices(C.c,c));

% only radius
C = capsule.generateRandom('Radius',r);
assert(withinTol(C.r,r));

% dimension and center
C = capsule.generateRandom('Dimension',n,'Center',c);
assert(dim(C) == n && compareMatrices(C.c,c));

% dimension and radius
C = capsule.generateRandom('Dimension',n,'Radius',r);
assert(dim(C) == n && withinTol(C.r,r));

% center and radius
C = capsule.generateRandom('Center',c,'Radius',r);
assert(compareMatrices(C.c,c) && withinTol(C.r,r));

% dimension, center, and radius
C = capsule.generateRandom('Dimension',n,'Center',c,'Radius',r);
assert(dim(C) == n && compareMatrices(C.c,c) && withinTol(C.r,r));


% unify results
res = true;

% ------------------------------ END OF CODE ------------------------------
