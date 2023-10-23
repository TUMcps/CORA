function res = test_ellipsoid_isequal
% test_ellipsoid_isequal - unit test function of isequal
%
% Syntax:
%    res = test_ellipsoid_isequal
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

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate ellipsoid
E1 = ellipsoid([1 0;0 2],[0; 1]);
E2 = ellipsoid([1 0;0 2],[0; 0]);
E3 = E1;

% compare results
res = ~isequal(E1,E2) && isequal(E1,E3);

% ------------------------------ END OF CODE ------------------------------
