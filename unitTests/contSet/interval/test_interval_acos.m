function res = test_interval_acos
% test_interval_acos - unit test function of sine for intervals,
%    overloaded 'acos()' function for intervals
%
% Syntax:
%    res = test_interval_acos
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

% Check special values:
% x     acos(x)
% -1    pi
% 0     pi/2
% 1     0
I = interval([-1,0],[1,0]);
I_acos = acos(I);
I_true = interval([0,pi/2],[pi,pi/2]);
assert(isequal(I_acos,I_true,tol));

% check out of bounds
I = interval(-2,0);
assertThrowsAs(@acos,'CORA:outOfDomain',I);

I = interval(0.5,1.1);
assertThrowsAs(@acos,'CORA:outOfDomain',I);

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
assertThrowsAs(@acos,'CORA:outOfDomain',I);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
