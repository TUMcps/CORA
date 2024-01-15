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

res = true(0);

% empty case
n = 2;
I = interval.empty(n);
res(end+1,1) = dim(I) == n;

% bounded vector
I = interval([-3; -2; -5],[4; 2; 1]);
n = dim(I);
n_true = 3;
res(end+1,1) = n == n_true;

% bounded matrix
I = interval([-3 0; -2 -1; -5 -2],[4 1; 2 3; 1 1]);
n = dim(I);
n_true = [3, 2];
res(end+1,1) = all(n == n_true);

% unbounded vector
I = interval([-Inf; -2; 3],[4; Inf; Inf]);
n = dim(I);
n_true = 3;
res(end+1,1) = n == n_true;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
