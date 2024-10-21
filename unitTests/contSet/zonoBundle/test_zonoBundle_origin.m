function res = test_zonoBundle_origin
% test_zonoBundle_origin - unit test function of origin
%
% Syntax:
%    res = test_zonoBundle_origin
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
zB = zonoBundle.origin(1);
zB_true = zonoBundle({zonotope(0)});
assert(isequal(zB,zB_true));
assert(contains(zB,0));

% 2D
zB = zonoBundle.origin(2);
zB_true = zonoBundle({zonotope(zeros(2,1))});
assert(isequal(zB,zB_true));
assert(contains(zB,zeros(2,1)));

% wrong calls
assertThrowsAs(@zonoBundle.origin,'CORA:wrongValue',0);
assertThrowsAs(@zonoBundle.origin,'CORA:wrongValue',-1);
assertThrowsAs(@zonoBundle.origin,'CORA:wrongValue',0.5);
assertThrowsAs(@zonoBundle.origin,'CORA:wrongValue',[1,2]);
assertThrowsAs(@zonoBundle.origin,'CORA:wrongValue','text');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
