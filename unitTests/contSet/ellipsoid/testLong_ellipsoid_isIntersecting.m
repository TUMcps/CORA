function res = testLong_ellipsoid_isIntersecting
% testLong_ellipsoid_isIntersecting - unit test function of
%    isIntersecting
%
% Syntax:
%    res = testLong_ellipsoid_isIntersecting
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
% Written:       18-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% point
[~,res] = evalc('testLong_ellipsoid_contains');

% all that implement distance
[~,res_] = evalc('testLong_ellipsoid_distance');
res = res && res_;

% mixed
res = res && testLong_component_ellipsoid_isIntersectingMixed;

% ------------------------------ END OF CODE ------------------------------
