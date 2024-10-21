function res = test_interval_min
% test_interval_min - unit test function of min
%
% Syntax:
%    res = test_interval_min
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

% empty interval, empty intervals
I = interval.empty(2);
I_min = min(I, I);
assert(representsa_(I_min,'emptySet',eps));

% empty interval, empty numeric
I = interval.empty(2);
I_min = min(I, zeros(2,0));
assert(representsa_(I_min,'emptySet',eps));

% bounded interval, numeric
I1 = interval([-2;-1],[2;1]);
I_min = min(I1, 0);
assert(isequal(I_min, interval([-2;-1], [0;0])));

% bounded interval, Inf
I1 = interval([-2;-1],[2;1]);
I_min = min(I1, Inf);
assert(isequal(I_min, I1));

% bounded interval, bounded interval
I1 = interval([-2;-1],[2;1]);
I2 = interval([-1;1], [1;3]);
I_min = min(I1, I2);
assert(isequal(I_min, interval([-2;-1], [1;1])));

% bounded interval, bounded interval (matrix)
I1 = interval([-2 0;-1 1],[2 2;1 4]);
I2 = interval([-1 -1;1 2],[1 -1;3 6]);
I_min = min(I1, I2);
assert(isequal(I_min, interval([-2 -1;-1 1],[1 -1;1 4])));

% bounded interval, unbounded interval
I1 = interval([-2;-1],[2;1]);
I2 = interval([-Inf;2],[-1;Inf]);
I_min = min(I1, I2);
assert(isequal(I_min, interval([-Inf;-1], [-1;1])));

% bounded interval, zonotope
I1 = interval([-2;-1],[2;1]);
Z = zonotope([1;-1], [2;1]);
I_min = min(I1, Z);
assert(isequal(I_min, interval([-2;-2], [2;0])));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);

v = 5;
I_min= min(I,v);
assert(max(I_min.inf >= v,[],'all'));

I_min = min(I,[],1);
assert(isequal(size(I_min),[1,2,2,3]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
