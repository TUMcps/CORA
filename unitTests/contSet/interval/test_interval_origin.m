function res = test_interval_origin
% test_interval_origin - unit test function of origin
%
% Syntax:
%    res = test_interval_origin
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
I = interval.origin(1);
I_true = interval(0);
assert(isequal(I,I_true));
assert(contains(I,0));

% 2D
I = interval.origin(2);
I_true = interval(zeros(2,1));
assert(isequal(I,I_true));
assert(contains(I,zeros(2,1)));

% wrong calls
assertThrowsAs(@interval.origin,'CORA:wrongValue',0);
assertThrowsAs(@interval.origin,'CORA:wrongValue',-1);
assertThrowsAs(@interval.origin,'CORA:wrongValue',0.5);
assertThrowsAs(@interval.origin,'CORA:wrongValue',[1,2]);
assertThrowsAs(@interval.origin,'CORA:wrongValue','text');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
