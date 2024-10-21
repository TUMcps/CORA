function res = test_interval_cat
% test_interval_cat - unit test function of cat
%
% Syntax:
%    res = test_interval_cat
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

% Authors:       Tobias Ladner
% Written:       18-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simple cases
I1 = interval([-0.5;0.3]);
I2 = interval([2;3]);
I = cat(2,I1,I2);
assert(isequal(dim(I),[2,2]));

I = cat(1,I1,I2);
assert(isequal(dim(I),4));

% n-d case
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I1 = interval(lb,ub);
I2 = interval(lb,ub);

I = cat(5,I1,I2,I1);
assert(isequal(dim(I),[2 2 2 3 3]));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
