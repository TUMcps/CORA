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

res = true(0);

% instantiate capsule
C = capsule([1; 0], [1; -1], 2);

% linear map: 90 degrees = pi/2
angle = pi/2;
R=[cos(angle) -sin(angle); sin(angle) cos(angle)];
C_mtimes = R * C;

% true solution
C_mtimes_true = capsule([0; 1], [1; 1], 2);

% compare results
res(end+1,1) = compareMatrices(C_mtimes.c,C_mtimes_true.c);
res(end+1,1) = compareMatrices(C_mtimes.g,C_mtimes_true.g);
res(end+1,1) = withinTol(C_mtimes.r,C_mtimes_true.r);

% empty set
C = capsule.empty(2);
res(end+1,1) = representsa_(R * C,'emptySet',eps);

% add results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
