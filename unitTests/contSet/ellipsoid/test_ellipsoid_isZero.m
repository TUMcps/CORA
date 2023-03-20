function res = test_ellipsoid_isZero
% test_ellipsoid_isZero - unit test function of isZero
%
% Syntax:  
%    res = test_ellipsoid_isZero
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

% Author:       Mark Wetzlinger
% Written:      17-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
res(1) = ~isZero(ellipsoid());

% only origin
E = ellipsoid(zeros(3),zeros(3,1));
res(2) = isZero(E);

% shifted center
E = ellipsoid(zeros(3),0.01*ones(3,1));
res(3) = ~isZero(E);

% shifted center, contains origin within tolerance
E = ellipsoid(0.01*eye(3),0.01*ones(3,1));
tol = 0.15;
res(4) = isZero(E,tol);

% combine results
res = all(res);

%------------- END OF CODE --------------
