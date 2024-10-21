function res = test_interval_dim
% test_interval_dim - unit test function of dim
%
% Syntax:
%    res = test_interval_dim
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
% Written:       27-July-2021
% Last update:   03-December-2023 (MW, add unbounded and matrix case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
n = 2;
I = interval.empty(n);
assert(dim(I) == n);

% bounded vector
I = interval([-3; -2; -5],[4; 2; 1]);
n = dim(I);
n_true = 3;
assert(n == n_true);

% bounded matrix
I = interval([-3 0; -2 -1; -5 -2],[4 1; 2 3; 1 1]);
n = dim(I);
n_true = [3, 2];
assert(all(n == n_true));

% unbounded vector
I = interval([-Inf; -2; 3],[4; Inf; Inf]);
n = dim(I);
n_true = 3;
assert(n == n_true);

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(isequal(dim(I),[2 2 2 3]))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
