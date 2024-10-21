function res = test_levelSet_dim
% test_levelSet_dim - unit test function of dim
%
% Syntax:
%    res = test_levelSet_dim
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
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init symbolic variables
syms x y

% empty level set
n = 2;
ls = levelSet.empty(n);
assert(dim(ls) == n);

% single equation
eq = x^2 + y^2 - 4;
ls = levelSet(eq,[x;y],'==');
assert(dim(ls) == 2);

% multiple equations with same variables
eq1 = x - y;
eq2 = x^2*y;
ls = levelSet([eq1;eq2],[x;y],{'<=','<'});
assert(dim(ls) == 2);

% equation with unused variable
eq = x;
ls = levelSet(eq,[x;y],'==');
assert(dim(ls) == 2);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
