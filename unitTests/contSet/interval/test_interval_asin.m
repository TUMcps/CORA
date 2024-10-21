function res = test_interval_asin
% test_interval_asin - unit test function of sine for intervals,
%    overloaded 'asin()' function for intervals
%
% Syntax:
%    res = test_interval_asin
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
% See also: mtimes

% Authors:       Mark Wetzlinger
% Written:       08-August-2020
% Last update:   03-December-2023 (MW, add out of bounds cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-9;

% Check special values
% x     asin(x)
% -1    -pi/2
% 0     0
% 1     pi/2

I = interval([-1,0],[1,0]);
I_asin = asin(I);
I_true = interval([-pi/2,0],[pi/2,0]);
assert(isequal(I_asin,I_true,tol));

% check out of bounds
I = interval(-2,0);
assertThrowsAs(@asin,'CORA:outOfDomain',I);

I = interval(0.5,1.1);
assertThrowsAs(@asin,'CORA:outOfDomain',I);

% n-d arrays
lb = [];
lb(:,:,1,1) = [1 2; 3 5];
lb(:,:,1,2) = [0 -1; -2 3];
lb(:,:,1,3) = [1 1; -1 0];
lb(:,:,2,1) = [-3 2; 0 1];
ub = [];
ub(:,:,1,1) = [1.5 4; 4 10];
ub(:,:,1,2) = [1 2; 0 4];
ub(:,:,1,3) = [2 3; -0.5 2];
ub(:,:,2,1) = [-1 3; 0 2];
I = interval(lb,ub);
assertThrowsAs(@asin,'CORA:outOfDomain',I);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
