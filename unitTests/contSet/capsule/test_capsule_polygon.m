function res = test_capsule_polygon
% test_capsule_polygon - unit test function of polygon
%
% Syntax:
%    res = test_capsule_polygon
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

% Authors:       Mark Wetzlinger
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init 2D capsule
C = capsule([1;-1],[4;2],1);

% compute polygon (inner-approximation)
p = polygon(C);

% check whether all points are contained in the capsule
res = all(contains(C,p));

% ------------------------------ END OF CODE ------------------------------
