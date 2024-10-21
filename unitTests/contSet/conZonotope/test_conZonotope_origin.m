function res = test_conZonotope_origin
% test_conZonotope_origin - unit test function of origin
%
% Syntax:
%    res = test_conZonotope_origin
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
cZ = conZonotope.origin(1);
cZ_true = conZonotope(0);
assert(isequal(cZ,cZ_true));
assert(contains(cZ,0));

% 2D
cZ = conZonotope.origin(2);
cZ_true = conZonotope(zeros(2,1));
assert(isequal(cZ,cZ_true));
assert(contains(cZ,zeros(2,1)));

% wrong calls
assertThrowsAs(@conZonotope.origin,'CORA:wrongValue',0);
assertThrowsAs(@conZonotope.origin,'CORA:wrongValue',-1);
assertThrowsAs(@conZonotope.origin,'CORA:wrongValue',0.5);
assertThrowsAs(@conZonotope.origin,'CORA:wrongValue',[1,2]);
assertThrowsAs(@conZonotope.origin,'CORA:wrongValue','text');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
