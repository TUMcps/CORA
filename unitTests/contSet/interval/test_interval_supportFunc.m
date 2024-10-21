function res = test_interval_supportFunc
% test_interval_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = test_interval_supportFunc
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
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
I = interval.empty(1);
dir = 1;
[val,x] = supportFunc(I,dir,'upper');
assert(val == -Inf && isempty(x));
[val,x] = supportFunc(I,dir,'lower');
assert(val == Inf && isempty(x));

% 2D set
I = interval([-2;-1],[3;4]);
dir = [2;-1];
[val,x] = supportFunc(I,dir,'upper');
assert(val == 7 && all(x == [3;-1]));
[val,x] = supportFunc(I,dir,'lower');
assert(val == -8 && all(x == [-2;4]));
[val,x] = supportFunc(I,dir,'range');
assert(isequal(val,interval(-8,7)) && all(all(x == [-2 3; 4 -1])));

% unbounded set
I = interval([2;-Inf;-1],[Inf;4;2]);
dir = [-2;1;-1];
[val,x] = supportFunc(I,dir);
assert(val == 1 && all(x == [2;4;-1]));
dir = [1;-2;1];
[~,x] = supportFunc(I,dir);
assert(all(x == [Inf;-Inf;2]));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
