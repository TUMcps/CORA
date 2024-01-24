function res = test_levelSet_empty
% test_levelSet_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_levelSet_empty
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
ls = levelSet.empty(n);
res(end+1,1) = representsa(ls,'emptySet') && dim(ls) == 1;

% 5D
n = 5;
ls = levelSet.empty(n);
res(end+1,1) = representsa(ls,'emptySet') && dim(ls) == 5;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
