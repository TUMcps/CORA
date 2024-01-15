function res = test_interval_length
% test_interval_length - unit_test_function of length
%
% Syntax:
%    res = test_interval_length
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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       19-January-2016
% Last update:   04-December-2023 (MW, add unbounded and matrix cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty
I = interval.empty(2);
res(end+1,1) = length(I) == 0;

% bounded
I = interval([-5, -4, -3, 0, 0, 5], [-2, 0, 2, 0, 5, 8]);
res(end+1,1) = length(I) == 6;
I = interval([-5; -4; -3; 0; 0; 5], [-2; 0; 2; 0; 5; 8]);
res(end+1,1) = length(I) == 6;

% bounded, matrix
I = interval([1 2 3; -2 1 -1],[3 4 6; -1 2 0]);
res(end+1,1) = length(I) == 3;

% unbounded
I = interval([-Inf;2],[1;Inf]);
res(end+1,1) = length(I) == 2;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
