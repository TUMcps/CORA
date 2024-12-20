function res = test_zonoBundle_empty
% test_zonoBundle_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_zonoBundle_empty
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

% Authors:       Tobias Ladner
% Written:       16-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D
n = 1;
zB = zonoBundle.empty(n);
assert(representsa(zB,'emptySet') && dim(zB) == 1);

% 5D
n = 5;
zB = zonoBundle.empty(n);
assert(representsa(zB,'emptySet') && dim(zB) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
