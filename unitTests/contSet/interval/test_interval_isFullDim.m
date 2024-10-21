function res = test_interval_isFullDim
% test_interval_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_interval_isFullDim
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
% Written:       27-July-2021
% Last update:   04-December-2023 (MW, add unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
I = interval.empty(2);
assert(~isFullDim(I));

% bounded, full-dimensional
I = interval([-2; -4; -7; -1; -2],[4; 2; 6; 4; 8]);
assert(isFullDim(I));

% bounded, degenerate
I = interval([-2; -4; 0; -1; -2],[4; 2; 0; 4; 8]);
assert(~isFullDim(I));

% unbounded, full-dimensional
I = interval([-Inf;-2],[1;1]);
assert(isFullDim(I));

% unbounded, degenerate
I = interval([-Inf;0],[1;0]);
assert(~isFullDim(I));

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(~isFullDim(I))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
