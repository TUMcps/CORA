function res = test_capsule_mtimes
% test_capsule_mtimes - unit test function of mtimes
%
% Syntax:
%    res = test_capsule_mtimes
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

% instantiate capsule
C = capsule([1; 0], [1; -1], 2);

% linear map: 90 degrees = pi/2
angle = pi/2;
R=[cos(angle) -sin(angle); sin(angle) cos(angle)];
C_mtimes = R * C;

% true solution
C_mtimes_true = capsule([0; 1], [1; 1], 2);

% compare results
assert(compareMatrices(C_mtimes.c,C_mtimes_true.c));
assert(compareMatrices(C_mtimes.g,C_mtimes_true.g));
assert(withinTol(C_mtimes.r,C_mtimes_true.r));

% empty set
C = capsule.empty(2);
assert(representsa_(R * C,'emptySet',eps));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
