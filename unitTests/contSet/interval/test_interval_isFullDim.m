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

res = true(0);

% empty case
I = interval.empty(2);
res(end+1,1) = ~isFullDim(I);

% bounded, full-dimensional
I = interval([-2; -4; -7; -1; -2],[4; 2; 6; 4; 8]);
res(end+1,1) = isFullDim(I);

% bounded, degenerate
I = interval([-2; -4; 0; -1; -2],[4; 2; 0; 4; 8]);
res(end+1,1) = ~isFullDim(I);

% unbounded, full-dimensional
I = interval([-Inf;-2],[1;1]);
res(end+1,1) = isFullDim(I);

% unbounded, degenerate
I = interval([-Inf;0],[1;0]);
res(end+1,1) = ~isFullDim(I);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
