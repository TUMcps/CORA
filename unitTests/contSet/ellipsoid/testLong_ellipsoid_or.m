function res = testLong_ellipsoid_or
% testLong_ellipsoid_or - unit test function of or
%
% Syntax:
%    res = testLong_ellipsoid_or
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
% Written:       14-October-2019
% Last update:   16-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% double
res = res && testLong_component_ellipsoid_orDouble();

% ellipsoid
res = res && testLong_component_ellipsoid_orEllipsoid();

% ------------------------------ END OF CODE ------------------------------
