function res = test_ellipsoid_copy
% test_ellipsoid_copy - unit test function of copy
%
% Syntax:
%    res = test_ellipsoid_copy
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
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D ellipsoid
Q = [2 0; 0 1]; q = [1;-1];
E = ellipsoid(Q,q);
E_copy = copy(E);
assert(isequal(E,E_copy));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
