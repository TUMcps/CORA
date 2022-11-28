function res = testLongDuration_ellipsoid_minkDiff
% testLongDuration_ellipsoid_minkDiff - unit test function of minkDiff
%
% Syntax:  
%    res = testLongDuration_ellipsoid_minkDiff
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
res = res && testLongDuration_component_ellipsoid_minkDiffEllipsoid;

%------------- END OF CODE --------------
