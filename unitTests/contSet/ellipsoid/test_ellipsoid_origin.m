function res = test_ellipsoid_origin
% test_ellipsoid_origin - unit test function of origin
%
% Syntax:
%    res = test_ellipsoid_origin
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
% Written:       21-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D
E = ellipsoid.origin(1);
E_true = ellipsoid(0,0);
assert(isequal(E,E_true));
assert(contains(E,0));

% 2D
E = ellipsoid.origin(2);
E_true = ellipsoid(zeros(2,2),zeros(2,1));
assert(isequal(E,E_true));
assert(contains(E,zeros(2,1)));

% wrong calls
assertThrowsAs(@ellipsoid.origin,'CORA:wrongValue',0);
assertThrowsAs(@ellipsoid.origin,'CORA:wrongValue',-1);
assertThrowsAs(@ellipsoid.origin,'CORA:wrongValue',0.5);
assertThrowsAs(@ellipsoid.origin,'CORA:wrongValue',[1,2]);
assertThrowsAs(@ellipsoid.origin,'CORA:wrongValue','text');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
