function res = test_levelSet_Inf
% test_levelSet_Inf - unit test function of R^n instantiation
%
% Syntax:
%    res = test_levelSet_Inf
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
ls = levelSet.Inf(n);
res(end+1,1) = representsa(ls,'fullspace') && dim(ls) == 1;

% 5D
n = 5;
ls = levelSet.Inf(n);
res(end+1,1) = representsa(ls,'fullspace') && dim(ls) == 5;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
