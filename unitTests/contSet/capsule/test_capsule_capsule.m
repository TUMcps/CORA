function res = test_capsule_capsule
% test_capsule_capsule - unit test function of capsule (constructor)
%
% Syntax:
%    res = test_capsule_capsule
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

% Authors:       Mark Wetzlinger
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% 2D center, generator, and radius
c = [2; 0];
g = [1; -1];
r = 0.5;

% admissible initializations
% only center
C = capsule(c);
assert(compareMatrices(C.c,c) && compareMatrices(C.g,zeros(2,1)) && withinTol(C.r,0));

% center and generator
C = capsule(c,g);
assert(compareMatrices(C.c,c) && compareMatrices(C.g,g));

% center, generator, and radius
C = capsule(c,g,r);
assert(compareMatrices(C.c,c) && compareMatrices(C.g,g) && withinTol(C.r,r));

% wrong initializations
cn_1 = [3; 2; -1];
gn_1 = [-2; 1; 1];
rneg = -0.2;
rvec = [1; 4];

% mismatch between center and generator
assertThrowsAs(@capsule,'CORA:wrongInputInConstructor',c,gn_1);
assertThrowsAs(@capsule,'CORA:wrongInputInConstructor',cn_1,g);
    
% negative radius
assertThrowsAs(@capsule,'CORA:wrongValue',c,g,rneg);

% radius as a vector
assertThrowsAs(@capsule,'CORA:wrongValue',c,g,rvec);

% too many input arguments
assertThrowsAs(@capsule,'CORA:numInputArgsConstructor',c,g,r,r);

% ------------------------------ END OF CODE ------------------------------
