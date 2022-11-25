function res = testLongDuration_ellipsoid_and
% test_ellipsoid_and - unit test function of and
%
% Syntax:  
%    res = test_ellipsoid_and
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
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

if res
    disp('test_ellipsoid_and successful');
else
    disp('test_ellipsoid_and failed');
end
%------------- END OF CODE --------------
