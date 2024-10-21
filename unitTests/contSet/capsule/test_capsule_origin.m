function res = test_capsule_origin
% test_capsule_origin - unit test function of origin
%
% Syntax:
%    res = test_capsule_origin
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
C = capsule.origin(1);
C_true = capsule(0);
assert(isequal(C,C_true));
assert(contains(C,0));

% 2D
C = capsule.origin(2);
C_true = capsule(zeros(2,1));
assert(isequal(C,C_true));
assert(contains(C,zeros(2,1)));

% wrong calls
assertThrowsAs(@capsule.origin,'CORA:wrongValue',0);
assertThrowsAs(@capsule.origin,'CORA:wrongValue',-1);
assertThrowsAs(@capsule.origin,'CORA:wrongValue',0.5);
assertThrowsAs(@capsule.origin,'CORA:wrongValue',[1,2]);
assertThrowsAs(@capsule.origin,'CORA:wrongValue','text');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
