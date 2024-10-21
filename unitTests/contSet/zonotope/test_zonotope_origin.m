function res = test_zonotope_origin
% test_zonotope_origin - unit test function of origin
%
% Syntax:
%    res = test_zonotope_origin
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
Z = zonotope.origin(1);
Z_true = zonotope(0);
assert(isequal(Z,Z_true));
assert(contains(Z,0));

% 2D
Z = zonotope.origin(2);
Z_true = zonotope(zeros(2,1));
assert(isequal(Z,Z_true));
assert(contains(Z,zeros(2,1)));

% wrong calls
assertThrowsAs(@zonotope.origin,'CORA:wrongValue',0);
assertThrowsAs(@zonotope.origin,'CORA:wrongValue',-1);
assertThrowsAs(@zonotope.origin,'CORA:wrongValue',0.5);
assertThrowsAs(@zonotope.origin,'CORA:wrongValue',[1,2]);
assertThrowsAs(@zonotope.origin,'CORA:wrongValue','text');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
