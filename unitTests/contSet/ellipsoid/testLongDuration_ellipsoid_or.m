function res = testLongDuration_ellipsoid_or
% testLongDuration_ellipsoid_or - unit test function of or
%
% Syntax:  
%    res = testLongDuration_ellipsoid_or
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
% Written:      14-October-2019
% Last update:  16-March-2021
% Last revision:---

%------------- BEGIN CODE --------------
res = true;

% double
res = res && testLongDuration_component_ellipsoid_orDouble();

% ellipsoid
res = res && testLongDuration_component_ellipsoid_orEllipsoid();


if res
    disp('testLongDuration_ellipsoid_or successful');
else
    disp('testLongDuration_ellipsoid_or failed');
end
%------------- END OF CODE --------------