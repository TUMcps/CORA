function res = test_capsule_plus
% test_capsule_plus - unit test function of plus
%
% Syntax:
%    res = test_capsule_plus
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
% Written:       28-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D capsule
c = [1; 0]; g = [-4; 3]; r = 1;
C = capsule(c,g,r);
% empty capsule
C_empty = capsule.empty(2);

% Minkowski sum: vector
vect = [1; 1];
C_vect = C + vect;
C_vect_true = capsule([2; 1], [-4; 3], 1);

res(end+1,1) = compareMatrices(C_vect.c,C_vect_true.c);
res(end+1,1) = compareMatrices(C_vect.g,C_vect_true.g);
res(end+1,1) = withinTol(C_vect.r,C_vect_true.r);

% Minkowski sum: capsule
C_add = capsule([0; 1], [sqrt(2); sqrt(2)], 1);
C_caps = C + C_add;
C_caps_true = capsule([1; 1], [-4; 3], 4);

res(end+1,1) = compareMatrices(C_caps.c,C_caps_true.c);
res(end+1,1) = compareMatrices(C_caps.g,C_caps_true.g);
res(end+1,1) = withinTol(C_caps.r,C_caps_true.r);

% Minkowski sum: empty set
res(end+1,1) = representsa_(C + C_empty,'emptySet',eps);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
