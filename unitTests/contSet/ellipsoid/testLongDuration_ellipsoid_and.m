function res = testLongDuration_ellipsoid_and
% testLongDuration_ellipsoid_and - unit test function of and
%
% Syntax:  
%    res = testLongDuration_ellipsoid_and
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

% Author:       Victor Gassmann
% Written:      17-October-2019
% Last update:  17-March-2021
% Last revision:---

%------------- BEGIN CODE --------------
res = true;

% ellipsoid
res = res && testLongDuration_component_ellipsoid_andEllipsoid();

% halfspace
res = res && testLongDuration_component_ellipsoid_andHalfspace();

%------------- END OF CODE --------------
