function res = test_capsule_uminus
% test_capsule_uminus - unit test function of uminus
%
% Syntax:
%    res = test_capsule_uminus
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

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate capsule
C = capsule([1;4], [1;-2], 2);

% negate
nC = -C;
assert(all(nC.c == -C.c));
assert(all(nC.g == -C.g));
assert(all(nC.r == C.r));

% compare with -1 * C
assert(isequal(nC, -1*C));

% second example
C = capsule([2;-3], [5;-6], 1);
nC = -C;
assert(isequal(nC, -1*C));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
