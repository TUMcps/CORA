function res = test_ellipsoid_isemptyobject
% test_ellipsoid_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_ellipsoid_isemptyobject
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty ellipsoid
E = ellipsoid.empty(2);
assert(isemptyobject(E));

% 2D ellipsoid
Q = eye(2); q = [-1;1];
E = ellipsoid(Q,q);
assert(~isemptyobject(E));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
