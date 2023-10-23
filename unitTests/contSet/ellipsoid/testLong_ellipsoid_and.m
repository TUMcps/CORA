function res = testLong_ellipsoid_and
% testLong_ellipsoid_and - unit test function of and
%
% Syntax:
%    res = testLong_ellipsoid_and
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

% Authors:       Victor Gassmann
% Written:       17-October-2019
% Last update:   17-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% ellipsoid
res = res && testLong_component_ellipsoid_andEllipsoid();

% halfspace
res = res && testLong_component_ellipsoid_andHalfspace();

% ------------------------------ END OF CODE ------------------------------
