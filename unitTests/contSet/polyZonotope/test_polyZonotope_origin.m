function res = test_polyZonotope_origin
% test_polyZonotope_origin - unit test function of origin
%
% Syntax:
%    res = test_polyZonotope_origin
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
% Written:       21-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D
pZ = polyZonotope.origin(1);
pZ_true = polyZonotope(0);
assert(isequal(pZ,pZ_true));
assert(contains(pZ,0));

% 2D
pZ = polyZonotope.origin(2);
pZ_true = polyZonotope(zeros(2,1));
assert(isequal(pZ,pZ_true));
assert(contains(pZ,zeros(2,1)));

% wrong calls
assertThrowsAs(@polyZonotope.origin,'CORA:wrongValue',0);
assertThrowsAs(@polyZonotope.origin,'CORA:wrongValue',-1);
assertThrowsAs(@polyZonotope.origin,'CORA:wrongValue',0.5);
assertThrowsAs(@polyZonotope.origin,'CORA:wrongValue',[1,2]);
assertThrowsAs(@polyZonotope.origin,'CORA:wrongValue','text');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
