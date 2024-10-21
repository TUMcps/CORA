function res = test_interval_empty
% test_interval_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_interval_empty
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
I = interval.empty(n);
assert(representsa(I,'emptySet') && dim(I) == 1);

% 5D
n = 5;
I = interval.empty(n);
assert(representsa(I,'emptySet') && dim(I) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
