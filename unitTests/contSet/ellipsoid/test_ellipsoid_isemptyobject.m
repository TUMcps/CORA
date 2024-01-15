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

res = true(0);

% empty ellipsoid
E = ellipsoid.empty(2);
res(end+1,1) = isemptyobject(E);

% 2D ellipsoid
Q = eye(2); q = [-1;1];
E = ellipsoid(Q,q);
res(end+1,1) = ~isemptyobject(E);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
