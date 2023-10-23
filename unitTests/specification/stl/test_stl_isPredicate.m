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
res = [];
res(end+1,1) = isPredicate(p1);
res(end+1,1) = isPredicate(p2);
res(end+1,1) = isPredicate(p3);
res(end+1,1) = isPredicate(p4);

res(end+1,1) = ~isPredicate(n1);
res(end+1,1) = ~isPredicate(n2);
res(end+1,1) = ~isPredicate(n3);
res(end+1,1) = ~isPredicate(n4);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
