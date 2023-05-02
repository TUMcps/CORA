function res = test_ellipsoid_display
% test_ellipsoid_display - unit test function of display
%
% Syntax:  
%    res = test_ellipsoid_display
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
% Written:      28-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% load ellipsoids
load cases.mat E_c

% empty set
E = ellipsoid()

% 2D ellipsoid
E = E_c{1}.E1

% degenerate ellipsoid
E = E_c{1}.Ed1

% all-zero shape matrix
E_c{1}.E0

%------------- END OF CODE --------------
