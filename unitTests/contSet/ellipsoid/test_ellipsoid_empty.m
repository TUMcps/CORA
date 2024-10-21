function res = test_ellipsoid_empty
% test_ellipsoid_empty - unit test function of empty instantiation
%
% Syntax:
%    res = test_ellipsoid_empty
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
E = ellipsoid.empty(n);
assert(representsa(E,'emptySet') && dim(E) == 1);

% 5D
n = 5;
E = ellipsoid.empty(n);
assert(representsa(E,'emptySet') && dim(E) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
