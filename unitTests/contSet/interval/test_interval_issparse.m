function res = test_interval_issparse
% test_interval_issparse - unit test function of issparse
%
% Syntax:
%    res = test_interval_issparse
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
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test both sparse
I = interval(speye(2), 2*speye(2));
assert(issparse(I));

% test either sparse
I = interval(speye(2), 2*eye(2));
assert(issparse(I));

I = interval(eye(2), 2*speye(2));
assert(issparse(I));

% test not sparse
I = interval([1;2], [3;4]);
assert(~issparse(I));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(~issparse(I))

res = true;


% ------------------------------ END OF CODE ------------------------------
