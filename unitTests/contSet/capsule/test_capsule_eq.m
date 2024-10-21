function res = test_capsule_eq
% test_capsule_eq - unit test function of '=='
%
% Syntax:
%    res = test_capsule_eq
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
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty capsule
C1 = capsule.empty(1);
C2 = capsule.empty(2);
assert(C1 == C1);
assert(~(C1 == C2));

% tolerance
tol = 1e-9;

% define properties
c = [2; 0; -1];
g = [0.2; -0.7; 0.4];
r = 1;
C = capsule(c,g,r);
C_tol = capsule(c,g,r+tol/2);
C_red = capsule(c(1:end-1),g(1:end-1),r);

% compare
assert(C == C);
assert(eq(C,C));
assert(~(C == C_tol));
assert(eq(C,C_tol,tol));
assert(~eq(C,C_red));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
