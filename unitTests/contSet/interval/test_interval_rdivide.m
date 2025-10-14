function res = test_interval_rdivide
% test_interval_rdivide - unit test function of ./ for intervals
%
% Syntax:
%    res = test_interval_rdivide
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

% Authors:       Tobias Ladner
% Written:       06-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I = interval(1,2);
assert(isequal(I./2,interval(0.5,1)))

I = interval([1 2; 3 4],[5 6; 7 8]);
assert(isequal(I./2,interval([0.5 1; 1.5 2],[2.5 3; 3.5 4])))

I = interval([1 2; 3 4],[5 6; 7 8]);
assert(isequal(I./[2 1],interval([0.5 2; 1.5 4],[2.5 6; 3.5 8])))

I = interval([1 2; 3 4],[5 6; 7 8]);
assert(isequal(I./interval([2 1]),interval([0.5 2; 1.5 4],[2.5 6; 3.5 8])))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
