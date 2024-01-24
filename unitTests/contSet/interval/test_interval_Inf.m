function res = test_interval_Inf
% test_interval_Inf - unit test function of R^n instantiation
%
% Syntax:
%    res = test_interval_Inf
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
I = interval.Inf(n);
res(end+1,1) = representsa(I,'fullspace') && dim(I) == 1;

% 5D
n = 5;
I = interval.Inf(n);
res(end+1,1) = representsa(I,'fullspace') && dim(I) == 5;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
