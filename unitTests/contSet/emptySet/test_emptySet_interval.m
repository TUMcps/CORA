function res = test_emptySet_interval
% test_emptySet_interval - unit test function of interval conversions
%
% Syntax:
%    res = test_emptySet_interval
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2d
O = emptySet(2);
I = interval(O);

assert(representsa(I,'emptySet'))
assert(dim(I) == 2)

% 3d
O = emptySet(3);
I = interval(O);

assert(representsa(I,'emptySet'))
assert(dim(I) == 3)

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
