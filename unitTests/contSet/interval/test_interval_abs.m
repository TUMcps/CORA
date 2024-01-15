function res = test_interval_abs
% test_interval_abs - unit test function of absolute value
%
% Syntax:
%    res = test_interval_abs
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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       05-January-2016
% Last update:   12-January-2016 (DG)
%                03-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-9;

% bounded and degenerate interval
I = interval([-5,-4,-3,0,0,5],[-2,0,2,0,5,8]);
Iabs = abs(I);
Iabs_true = interval([2,0,0,0,0,5],[5,4,3,0,5,8]);
res(end+1,1) = isequal(Iabs,Iabs_true,tol);

% unbounded interval
I = interval(2,Inf);
Iabs = abs(I);
res(end+1,1) = isequal(I,Iabs,tol);

I = interval([-Inf;-2],[-2;Inf]);
Iabs = abs(I);
Iabs_true = interval([2;0],[Inf;Inf]);
res(end+1,1) = isequal(Iabs,Iabs_true,tol);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
