function res = test_interval_randPoint
% test_interval_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_interval_randPoint
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       21-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
n = 2;
I = interval.empty(n);
p = randPoint(I);
assert(isempty(p) && isnumeric(p) && size(p,1) == n);

% 2D
lb = [-2; -1]; ub = [1; 0];
I = interval(lb,ub);
p = randPoint(I);
assert(contains(I,p));
p = randPoint(I,10);
assert(all(contains(I,p)));
p = randPoint(I,1,'extreme');
assert(contains(I,p));

% degenerate ---

% just center
lb = [-2; -1]; ub = [1; -1];
I = interval(lb,ub);
p = randPoint(I,10);
assert(all(contains(I,p)));
p = randPoint(I,1,'extreme');
assert(contains(I,p));

% only one dimension
lb = [-2; -1]; ub = [1; -1];
I = interval(lb,ub);
p = randPoint(I,10);
assert(all(contains(I,p)));
p = randPoint(I,1,'extreme');
assert(contains(I,p));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
