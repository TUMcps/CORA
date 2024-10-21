function res = test_capsule_empty
% test_capsule_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_capsule_empty
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
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D empty capsule
n = 2;
C = capsule.empty(n);
assert(dim(C) == n);
assert(representsa(C,'emptySet'));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
