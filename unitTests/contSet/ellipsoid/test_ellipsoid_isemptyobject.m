function res = test_ellipsoid_isemptyobject
% test_ellipsoid_isemptyobject - unit test function of isemptyobject
%
% Syntax:  
%    res = test_ellipsoid_isemptyobject
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
% Written:      03-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate ellipsoids
E1 = ellipsoid();

Q = eye(2);
q = [-1;1];
E2 = ellipsoid(Q,q);

% check results
res = isemptyobject(E1) && ~isemptyobject(E2);

%------------- END OF CODE --------------
