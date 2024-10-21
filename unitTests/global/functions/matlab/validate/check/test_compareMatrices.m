function res = test_compareMatrices
% test_compareMatrices - unit test function for comparsion of matrices
%
% Syntax:
%    res = test_compareMatrices()
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
% Written:       22-November-2022
% Last update:   08-May-2023 (TL, ordered)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init matrices
A = [1 2 3 4; 5 6 7 8; 9 10 11 12];
B = [1 2 3 4; 5 6 7 8; 9 10 11 12];
B_order = B(:,[4 2 3 1]);
B_eps = B_order + 10*eps;
B_minus = B_order(:,1:3);
B_plus = [B_order, [1; 1; 1]];
C = [9 10 11 12; 5 6 7 8; 1 2 3 4];
D = [];

% completely equal matrices
assert(compareMatrices(A,B))
assert(compareMatrices(A,B,eps,'subset'))

% compare with different order
assert(compareMatrices(A,B_order))
assert(compareMatrices(A,B_order,eps,'subset'))

% check with tolerance
assert(~compareMatrices(A,B_eps,eps))
assert(compareMatrices(A,B_eps,100*eps))

% check different sizes and subset
assert(~compareMatrices(A,B_plus))
assert(~compareMatrices(A,B_minus))

% check subset
assert(compareMatrices(B_minus,B,eps,'subset'))
assert(compareMatrices(A,B_plus,eps,'subset'))

% completely different matrices
assert(~compareMatrices(A,C))

% check ordered
assert(compareMatrices(A,B,eps,'equal',true))
assert(~compareMatrices(A(:,1:end-1),B,eps,'equal',true))
assert(compareMatrices(A,B,eps,'subset',true))
assert(compareMatrices(A(:,1:end-1),B,eps,'subset',true))
assert(~compareMatrices(A,B_order,eps,'equal',true))
assert(~compareMatrices(A,B_order,eps,'subset',true))


% check wrong input arguments

% matrices have to be nonempty
assertThrowsAs(@compareMatrices,'CORA:wrongValue',A,D);

% matrices have to be nonempty
assertThrowsAs(@compareMatrices,'CORA:wrongValue',D,B);

% eps has to be nonnegative
assertThrowsAs(@compareMatrices,'CORA:wrongValue',A,B,-1);

% eps has to be a scalar
assertThrowsAs(@compareMatrices,'CORA:wrongValue',A,B,[1 1]);

% flag has to be fourth input argument
assertThrowsAs(@compareMatrices,'CORA:wrongValue',A,B,'subset');

% flag has to be 'equal' or 'subset'
assertThrowsAs(@compareMatrices,'CORA:wrongValue',A,B,eps,'exact');


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
