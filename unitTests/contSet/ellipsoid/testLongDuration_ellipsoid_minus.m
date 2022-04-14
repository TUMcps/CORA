function res = testLongDuration_ellipsoid_minus
% testLongDuration_ellipsoid_minus - unit test function of minus
%
% Syntax:  
%    res = testLongDuration_ellipsoid_minus
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
% Written:      15-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;

% double
res = res && testLongDuration_component_ellipsoid_minusDouble;

% ellipsoid
res = res && testLongDuration_component_ellipsoid_minusEllipsoid;

if res
    disp('testLongDuration_ellipsoid_minus successful');
else
    disp('testLongDuration_ellipsoid_minus failed');
end
%------------- END OF CODE --------------
