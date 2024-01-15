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

res = true(0);

% empty ellipsoid
n = 2;
E = ellipsoid.empty(n);
c = center(E);
res(end+1,1) = isempty(c) && isnumeric(c) && size(c,1) == n;

% 2D ellipsoids
E1 = ellipsoid([1 0;0 2],[0; 1]);
E2 = ellipsoid([1 0;0 2]); % center at origin
res(end+1,1) = isequal(center(E1),[0;1]);
res(end+1,1) = isequal(center(E2),[0;0]);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
