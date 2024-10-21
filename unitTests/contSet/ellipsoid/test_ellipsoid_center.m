function res = test_ellipsoid_center
% test_ellipsoid_center - unit test function of center
%
% Syntax:
%    res = test_ellipsoid_center
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
% Written:       27-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty ellipsoid
n = 2;
E = ellipsoid.empty(n);
c = center(E);
assert(isempty(c) && isnumeric(c) && size(c,1) == n);

% 2D ellipsoids
E1 = ellipsoid([1 0;0 2],[0; 1]);
E2 = ellipsoid([1 0;0 2]); % center at origin
assert(isequal(center(E1),[0;1]));
assert(isequal(center(E2),[0;0]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
