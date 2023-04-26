function res = testLong_ellipsoid_minkDiff
% testLong_ellipsoid_minkDiff - unit test function of minkDiff
%
% Syntax:  
%    res = testLong_ellipsoid_minkDiff
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
% Written:      15-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% double
res = res && testLong_component_ellipsoid_minusDouble;

% ellipsoid
res = res && testLong_component_ellipsoid_minkDiffEllipsoid;

%------------- END OF CODE --------------
