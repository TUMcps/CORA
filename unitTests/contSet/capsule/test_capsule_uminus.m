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

resvec = true(0);

% instantiate capsule
C = capsule([1;4], [1;-2], 2);

% negate
nC = -C;
resvec(end+1) = all(nC.c == -C.c);
resvec(end+1) = all(nC.g == -C.g);
resvec(end+1) = all(nC.r == C.r);

% compare with -1 * C
resvec(end+1) = isequal(nC, -1*C);

% second example
C = capsule([2;-3], [5;-6], 1);
nC = -C;
resvec(end+1) = isequal(nC, -1*C);

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
