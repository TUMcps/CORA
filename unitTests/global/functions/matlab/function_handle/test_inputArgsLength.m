function res = test_inputArgsLength
% test_inputArgsLength - unit test function for automated read out of
%    number of input arguments to a function handle
%
% Syntax:
%    res = test_inputArgsLength()
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
% Written:       20-November-2022
% Last update:   19-January-2024 (MW, test reported bug)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% for comments:
%    n ... state dimension,
%    m ... input dimension,
%    r ... output dimension.

% empty function handle
f = @(x,u) [];
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); m = inputArgs(2);
assert(n == 0);
assert(m == 0);
assert(r == 0);

% 1. One input variable: x

% n = 1, r = 1
f = @(x) x(1)^2;
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1);
assert(length(inputArgs) == 1 && n == 1 && r == 1);

% n = 2, r = 1
f = @(x) x(1)^2 + x(2) - 4;
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1);
assert(length(inputArgs) == 1 && n == 2 && r == 1);

% n = 5, r = 3
f = @(x) [x(1)^2 + x(5)/5; 2*x(2) + x(3); x(4)^3];
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1);
assert(length(inputArgs) == 1 && n == 5 && r == 3);

% n = 5, r = 3 (but not all x are used)
f = @(x) [x(1)^2 + x(5)/3; x(2)*2; 1-x(4)^3];
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1);
assert(length(inputArgs) == 1 && n == 5 && r == 3);

% n = 4, r = 2 (matrix multiplication with indexing)
f = @(x) randn(2,4) * x(1:4);
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1);
assert(length(inputArgs) == 1 && n == 4 && r == 2);

% n = 5, r = 2 (matrix multiplication with zeros)
M = [0 1 0 0 0; 0 0 1 0 0];
f = @(x) M * x(1:5);
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1);
assert(length(inputArgs) == 1 && n == 5 && r == 2);

% n = 4, r = 2 (matrix multiplication without indexing)
f = @(x) randn(2,4) * x;
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1);
assert(length(inputArgs) == 1 && n == 4 && r == 2);

% n = 4, r = 2 (matrix addition without indexing)
f = @(x) randn(2,1) + x;
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1);
% This will most likely fail, because the function definition is ambiguous;
% should x be a vector? A scalar added component-wise?
% For now, we negate this test; if you see that it gives you 'false' here,
% then either you did something wrong, or you found a better way to deal
% with this problem
assert(length(inputArgs) == 1 && ~(n == 2 && r == 2));

% 2. Two input variables: x, u

% n = 1, m = 1, r = 1
f = @(x,u) x(1)^2 + u(1);
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); m = inputArgs(2);
assert(length(inputArgs) == 2 && n == 1 && m == 1 && r == 1);

% n = 6, m = 4, r = 2
f = @(x,u) [x(1)^2 + u(1); u(4)^2 - x(6)];
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); m = inputArgs(2);
assert(length(inputArgs) == 2 && n == 6 && m == 4 && r == 2);

% n = 4, m = 2, r = 2 (matrix multiplication with indexing)
f = @(x,u) randn(2,4) * x(1:4) + randn(2,2) * u(1:2);
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); m = inputArgs(2);
assert(length(inputArgs) == 2 && n == 4 && m == 2 && r == 2);

% n = 4, r = 2 (matrix multiplication without indexing)
f = @(x,u) randn(2,4) * x + randn(2,2) * u;
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); m = inputArgs(2);
assert(length(inputArgs) == 2 && n == 4 && m == 2 && r == 2);

% n = 5, m = 2, r = 2 (matrix multiplication with zeros)
M = [0 1 0 0 0; 0 0 1 0 0];
N = [1 0; 0 0];
f = @(x,u) M * x(1:5) + N * u(1:2);
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); m = inputArgs(2);
assert(length(inputArgs) == 2 && n == 5 && m == 2 && r == 2);

% n = 2, m = 0, r = 2
f = @(x,u) [x(1)^2 - x(2), x(2)];
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); m = inputArgs(2);
assert(length(inputArgs) == 2 && n == 2 && m == 0 && r == 2);

% 3. Three input variables: x, y, u

% n = 3, m1 = 2, m2 = 1, r = 2
f = @(x,y,u) [x(1)^2 + y(1) - x(2), x(2)*x(3) + u(1) - y(2)^2];
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); y = inputArgs(2); m = inputArgs(3);
assert(length(inputArgs) == 3 && n == 3 && y == 2 && m == 1 && r == 2);

% n = 2, m1 = 0, m2 = 1, r = 2
f = @(x,y,u) [x(1)^2 - x(2), x(2) + u(1)];
[inputArgs,r] = inputArgsLength(f);
n = inputArgs(1); y = inputArgs(2); m = inputArgs(3);
assert(length(inputArgs) == 3 && n == 2 && y == 0 && m == 1 && r == 2);

% n = 6, m = 2, o = 1
f = @(x,u,w) [x(4)*cos(x(3)) - x(5)*sin(x(3));
              x(4)*sin(x(3)) + x(5)*cos(x(3));
              x(6);
              x(5)*x(6) - sin(x(3)) + cos(x(3))*w(1);
              -x(5)*x(6) - cos(x(3)) + u(1) + u(2) - sin(x(3))*w(1);
              u(1) - u(2)];
inputArgs = inputArgsLength(f,3);
n = inputArgs(1); m = inputArgs(2); o = inputArgs(3);
assert(length(inputArgs) == 3 && n == 6 && m == 2 && o == 1);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
