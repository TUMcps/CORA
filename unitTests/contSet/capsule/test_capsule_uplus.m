function res = test_capsule_uplus
% test_capsule_uplus - unit test function of uminus
%
% Syntax:
%    res = test_capsule_uplus
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

% +
pC = +C;
assert(all(pC.c == C.c));
assert(all(pC.g == C.g));
assert(all(pC.r == C.r));

% compare with C
assert(isequal(pC, C));

% second example
C = capsule([2;-3], [5;-6], 1);
pC = +C;
assert(isequal(pC, +1*C));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
