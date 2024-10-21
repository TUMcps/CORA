function res = test_stl_isPredicate
% test_stl_isPredicate - unit test of stl isPredicate method
%
% Syntax:
%    res = test_stl_isPredicate
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Benedikt Seidl
% Written:       12-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

x = stl('x', 3);

% cases
p1 = x(1) < 5;
p2 = x(1) < 4 & x(2) > 7;
p3 = x(1) * x(2) + 4 < x(3) ^ 7;
p4 = stl('y');

n1 = finally(x(1) < 5, interval(1, 2));
n2 = ~(x(1) < 5);
n3 = x(1) < 5 & until(x(2) < 3, x(3) > 4, interval(1, 2));
n4 = stl(false);


% test
res = true;
assert(isPredicate(p1));
assert(isPredicate(p2));
assert(isPredicate(p3));
assert(isPredicate(p4));

assert(~isPredicate(n1));
assert(~isPredicate(n2));
assert(~isPredicate(n3));
assert(~isPredicate(n4));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
