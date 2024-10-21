function res = test_withinTol
% test_withinTol - unit test function for checking equality of values
%    within tolerance
%
% Syntax:
%    res = test_withinTol
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
% Written:       30-April-2023
% Last update:   03-December-2023 (MW, add Inf cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% scalar values
a = 3;
b = 3.1;
assert(withinTol(a,a));
assert(~withinTol(a,b));
assert(withinTol(a,b,b-a));

% scalar values (Inf)
a = -Inf;
b = Inf;
assert(withinTol(a,a));
assert(~withinTol(a,b));

% scalar vs. vector
a = 0;
b = [0.001 0.002 -0.001 0 0.003];
assert(all(size(withinTol(a,b)) == size(b)));
assert(all(withinTol(a,b,max(abs(b))-a)));

% scalar vs. vector (Inf)
a = -Inf;
b = [-Inf;-Inf];
assert(all(size(withinTol(a,b)) == size(b)));
assert(all(withinTol(a,b)));

% vector vs. vector
a = [3 4 5];
b = [4 5 6];
assert(all(size(withinTol(a,b)) == size(a)));
assert(~any(withinTol(a,b)));
assert(all(withinTol(a,b,1)));

% scalar vs. matrix
a = 0;
B = [0 0.001 0.002; -0.001 -0.004 0.003];
assert(all(size(withinTol(a,B)) == size(B)));
assert(all(all(withinTol(a,B,max(max(abs(B)))))));

% scalar vs. matrix (Inf)
a = Inf;
B = Inf(3,2);
assert(all(size(withinTol(a,B)) == size(B)));
assert(all(all(withinTol(a,B))));

% vector vs. matrix
a = [2 1];
B = [2 1.001; 1.999 0.999; 2.002 1.003];
assert(all(size(withinTol(a,B)) == size(B)));
assert(all(all(withinTol(a,B,max(max(abs(B-a)))))));

% vector vs. matrix (Inf)
a = [-Inf; Inf];
B = [-Inf -Inf; Inf Inf];
assert(all(size(withinTol(a,B)) == size(B)));
assert(all(all(withinTol(a,B))));

% matrix vs. matrix
A = [0 1; 1 0; 2 0];
B = [-0.001 1; 1.001 0.002; 1.999 -0.002];
assert(all(size(withinTol(A,B)) == size(A)));
assert(all(all(withinTol(A,B,max(max(abs(B-A)))))));

% wrong syntaxes
assertThrowsAs(@withinTol,'MATLAB:TooManyInputs',1,1,1,1);

% vector of different lengths
assertThrowsAs(@withinTol,'CORA:dimensionMismatch',[1 0],[1 0 0]);
assertThrowsAs(@withinTol,'CORA:dimensionMismatch',[1 0; 1 2],[1 0 0; 0 1 2]);

% tolerance not a scalar
assertThrowsAs(@withinTol,'CORA:wrongValue',[1 0],[1 0],[1 0]);

% tolerance not a nonnegative value
assertThrowsAs(@withinTol,'CORA:wrongValue',[1 0],[1 0],-1);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
