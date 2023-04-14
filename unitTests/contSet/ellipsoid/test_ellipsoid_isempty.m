function res = test_ellipsoid_isempty
% test_ellipsoid_isempty - unit test function of isempty
%
% Syntax:  
%    res = test_ellipsoid_isempty
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
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate ellipsoids
E1 = ellipsoid([1 0;0 2],[0; 1]);
% empty ellipsoids cannot be instantiated...

% compare results
res = ~isempty(E1) && isempty(ellipsoid());

%------------- END OF CODE --------------