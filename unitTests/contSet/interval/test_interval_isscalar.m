function res = test_interval_isscalar
% test_interval_isscalar - unit test function of isscalar
%
% Syntax:
%    res = test_interval_isscalar
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
% Written:       16-January-2016
% Last update:   04-December-2023 (MW, add unbounded and matrix cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty
I = interval.empty(1);
res(end+1,1) = ~isscalar(I);

% bounded, scalar
I = interval(-5, 2);
res(end+1,1) = isscalar(I);

% bounded, vector
I = interval([-5,-4,-3,0,0,5],[-2,0,2,0,5,8]);
res(end+1,1) = ~isscalar(I);

% bounded, matrix
I = interval([1 2; 3 4],[5 6; 7 8]);
res(end+1,1) = ~isscalar(I);

% unbounded, scalar
I = interval(-Inf,Inf);
res(end+1,1) = isscalar(I);

% unbounded, vector
I = interval([-Inf;1],[2;Inf]);
res(end+1,1) = ~isscalar(I);

% unbounded, matrix
I = interval([-Inf 2; 1 -Inf],[2 4; Inf 0]);
res(end+1,1) = ~isscalar(I);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
