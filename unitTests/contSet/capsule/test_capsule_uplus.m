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

resvec = true(0);

% instantiate capsule
C = capsule([1;4], [1;-2], 2);

% +
pC = +C;
resvec(end+1) = all(pC.c == C.c);
resvec(end+1) = all(pC.g == C.g);
resvec(end+1) = all(pC.r == C.r);

% compare with C
resvec(end+1) = isequal(pC, C);

% second example
C = capsule([2;-3], [5;-6], 1);
pC = +C;
resvec(end+1) = isequal(pC, +1*C);

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
