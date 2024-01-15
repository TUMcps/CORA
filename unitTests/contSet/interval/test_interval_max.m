function res = test_interval_max
% test_interval_max - unit test function of max
%
% Syntax:
%    res = test_interval_max
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

% Authors:       Tobias Ladner, Mark Wetzlinger
% Written:       16-December-2022
% Last update:   04-December-2023 (MW, add unbounded and matrix cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty interval, empty interval
I = interval.empty(2);
I_max = max(I, I);
res(end+1,1) = representsa(I_max,'emptySet');

% empty interval, empty numeric
I = interval.empty(2);
I_max = max(I, zeros(2,0));
res(end+1,1) = representsa(I_max,'emptySet');

% interval, numeric
I1 = interval([-2;-1],[2;1]);
I_max = max(I1, 0);
res(end+1,1) = isequal(I_max, interval([0;0], [2;1]));

% interval, interval
I1 = interval([-2;-1],[2;1]);
I2 = interval([-1;1], [1;3]);
I_max = max(I1, I2);
res(end+1,1) = isequal(I_max, interval([-1;1], [2;3]));

% interval, interval (matrix)
I1 = interval([-2 1;-1 3],[2 2;1 5]);
I2 = interval([-1 -4;1 2],[1 2;3 6]);
I_max = max(I1, I2);
res(end+1,1) = isequal(I_max, interval([-1 1;1 3],[2 2;3 6]));

% unbounded interval, bounded interval
I1 = interval([-Inf;2],[3;5]);
I2 = interval([-1;-3],[2;6]);
I_max = max(I1, I2);
res(end+1,1) = isequal(I_max, interval([-1;2],[3;6]));

% unbounded interval, bounded interval (unbounded result)
I1 = interval([-Inf;2],[Inf;5]);
I2 = interval([-1;-3],[Inf;6]);
I_max = max(I1, I2);
res(end+1,1) = isequal(I_max, interval([-1;2],[Inf;6]));

% interval, zonotope
I1 = interval([-2;-1],[2;1]);
Z = zonotope([1;-1], [2;1]);
I_max = max(I1, Z);
res(end+1,1) = isequal(I_max, interval([-1;-1], [3;1]));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
