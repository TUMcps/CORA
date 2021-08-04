function res = testLongDuration_ellipsoid_isIntersecting
% testLongDuration_ellipsoid_isIntersecting - unit test function of
%    isIntersecting
%
% Syntax:  
%    res = testLongDuration_ellipsoid_isIntersecting
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

% point
[~,res] = evalc('testLongDuration_ellipsoid_in');

% all that implement distance
[~,res_] = evalc('testLongDuration_ellipsoid_distance');
res = res && res_;

% mixed
res = res && testLongDuration_component_ellipsoid_isIntersectingMixed;

if res
    disp('testLongDuration_ellipsoid_isIntersecting successful');
else
    disp('testLongDuration_ellipsoid_isIntersecting failed');
end
%------------- END OF CODE --------------