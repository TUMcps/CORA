function res = testLong_ellipsoid_distance
% testLong_ellipsoid_distance - unit test function of distance
%
% Syntax:  
%    res = testLong_ellipsoid_distance
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
% Written:      18-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;

% ellipsoid
res = res && testLong_component_ellipsoid_distanceEllipsoid;

% conHyperplane
res = res && testLong_component_ellipsoid_distanceHyperplane;

% mptPolytope (halfspace is implicitly tested)
res = res && testLong_component_ellipsoid_distanceMptPolytope;

% double
res = res && testLong_component_ellipsoid_distancePoint;

%------------- END OF CODE --------------