function res = test_interval_diag
% test_interval_diag - unit test function of interval/diag
%
% Syntax:
%    res = test_interval_diag
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
% See also: interval/diag

% Authors:       Tobias Ladner
% Written:       18-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. Empty case
I = interval.empty(2);
assert(representsa(diag(I),'emptySet'));
assert(representsa(diag(I,0),'emptySet'));

    
% init random interval
lb = [-3; -2; -5];
ub = [4; 2; 1];
I = interval(lb,ub);

% test diagonal
D = diag(I);
assert(isequal(D,interval(diag(lb),diag(ub))));
assert(dim(I) == 3);

% test getter
assert(isequal(diag(D,0),I));
assert(isequal(diag(D,1),interval([0;0])));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assertThrowsAs(@diag,'CORA:wrongValue',I)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
