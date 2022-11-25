function res = testLongDuration_ellipsoid_distance
% testLongDuration_ellipsoid_distance - unit test function of distance
%
% Syntax:  
%    res = testLongDuration_ellipsoid_distance
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
% Written:      18-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;

% ellipsoid
res = res && testLongDuration_component_ellipsoid_distanceEllipsoid;

% conHyperplane
res = res && testLongDuration_component_ellipsoid_distanceHyperplane;

% mptPolytope (halfspace is implicitly tested)
res = res && testLongDuration_component_ellipsoid_distanceMptPolytope;

% double
res = res && testLongDuration_component_ellipsoid_distancePoint;

if res
    disp('testLongDuration_ellipsoid_distance successful');
else
    disp('testLongDuration_ellipsoid_distance failed');
end
%------------- END OF CODE --------------