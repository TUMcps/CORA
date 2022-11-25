function res = testLongDuration_ellipsoid_convHull
% testLongDuration_ellipsoid_convHull - unit test function of convHull
%
% Syntax:  
%    res = testLongDuration_ellipsoid_convHull
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

% simply is "or"
res = testLongDuration_ellipsoid_or;

if res
    disp('testLongDuration_ellipsoid_convHull successful');
else
    disp('testLongDuration_ellipsoid_convHull failed');
end
%------------- END OF CODE --------------