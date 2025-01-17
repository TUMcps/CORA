function res = test_interval_enclose
% test_interval_enclose - unit test function of the set enclosure
%
% Syntax:
%    res = test_interval_enclose
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

% Authors:       Maximilian Perschl
% Written:       18-November-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case: enclose
I_empty = interval.empty(3);
I_3D = interval(-rand(3,1),rand(3,1));
assert(isequal(enclose(I_3D,I_empty),I_3D));

% bounded, 1D
I1 = interval(-2,4);
% ...I2 contained in I1
I2 = interval(-1,3);
I_enclose = enclose(I1,I2);
assert(isequal(I_enclose,I1));
% ...I2 intersects I1
I2 = interval(3,6);
I_enclose = enclose(I1,I2);
I_true = interval(-2,6);
assert(isequal(I_enclose,I_true));
% ...I2 is disjoint from I1
I2 = interval(5,6);
I_enclose = enclose(I1,I2);
I_true = interval(-2,6);
assert(isequal(I_enclose,I_true));

% unbounded, 2D
I1 = interval([-2;-Inf],[2;1]);
I2 = interval([-1;3],[Inf;8]);
I_enclose = enclose(I1,I2);
I_true = interval([-2;-Inf],[Inf;8]);
assert(isequal(I_enclose,I_true));


% dimension mismatch
I_dim3 = interval(-rand(3,1),rand(3,1));
I_dim7 = interval(-rand(7,1),rand(7,1));
assertThrowsAs(@enclose,'CORA:dimensionMismatch',I_dim3,I_dim7);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
