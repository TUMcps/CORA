function res = test_interval_center
% test_interval_center - unit test function of center
%
% Syntax:
%    res = test_interval_center
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
% Written:       15-January-2016
% Last update:   03-December-2023 (MW, add unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true(0);

% 1. empty set
n = 2;
I = interval.empty(n);
c = center(I);
res(end+1,1) = isempty(c) && isnumeric(c) && all(size(c) == [2,0]);

% 2. bounded
I = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = center(I);
c_true = [-3.5,-2,-0.5,0,2.5,6.5];
res(end+1,1) = all(withinTol(c,c_true,tol));

% 3. unbounded
I = interval(-Inf,2);
c = center(I);
res(end+1,1) = isnan(c);
I = interval(2,Inf);
c = center(I);
res(end+1,1) = isnan(c);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
