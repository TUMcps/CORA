function res = test_emptySet_empty
% test_emptySet_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_emptySet_empty
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

res = true(0);

% 1D
n = 1;
O = emptySet.empty(n);
res(end+1,1) = representsa(O,'emptySet') && dim(O) == 1;

% 5D
n = 5;
O = emptySet.empty(n);
res(end+1,1) = representsa(O,'emptySet') && dim(O) == 5;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
