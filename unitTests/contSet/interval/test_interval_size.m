function res = test_interval_size
% test_interval_size - unit test function of size
%
% Syntax:
%    res = test_interval_size
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
% Written:       28-August-2019
% Last update:   05-December-2023 (MW, add empty/unbounded/matrix cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty
I = interval.empty(2);
s = size(I);
assert(all(size(s) == [1,2]) && all(s == [2,0]));

% vector, bounded
I = interval([-3;-2],[1;2]);
s = size(I);
assert(all(size(s) == [1,2]) && all(s == [2,1]));

% vector, unbounded
I = interval([-Inf,-2],[1,Inf]);
s = size(I);
assert(all(size(s) == [1,2]) && all(s == [1,2]));

% matrix
I = interval([-3 -2; 1 2; -3 -1],[4 5; 1 4; 0 1]);
s = size(I);
assert(all(size(s) == [1,2]) && all(s == [3,2]));

% n-d arrays
inf = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(inf,sup);
assert(isequal(size(I),[2,2,2,3]))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
