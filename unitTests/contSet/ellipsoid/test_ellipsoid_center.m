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

% Author:       Mark Wetzlinger
% Written:      27-September-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate ellipsoids
E1 = ellipsoid([1 0;0 2],[0; 1]);
E2 = ellipsoid([1 0;0 2]); % center at origin
E3 = ellipsoid();
res_e = isnumeric(center(E3)) && isempty(center(E3));
% compare results
res = isequal(center(E1),[0;1]) && isequal(center(E2),[0;0]) && res_e;

%------------- END OF CODE --------------