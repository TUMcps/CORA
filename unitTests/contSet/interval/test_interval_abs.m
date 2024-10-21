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

tol = 1e-9;

% bounded and degenerate interval
I = interval([-5,-4,-3,0,0,5],[-2,0,2,0,5,8]);
Iabs = abs(I);
Iabs_true = interval([2,0,0,0,0,5],[5,4,3,0,5,8]);
assert(isequal(Iabs,Iabs_true,tol));

% unbounded interval
I = interval(2,Inf);
Iabs = abs(I);
assert(isequal(I,Iabs,tol));

I = interval([-Inf;-2],[-2;Inf]);
Iabs = abs(I);
Iabs_true = interval([2;0],[Inf;Inf]);
assert(isequal(Iabs,Iabs_true,tol));

% n-d arrays
lb = [];
lb(:,:,1,1) = [1 2; 3 5];
lb(:,:,1,2) = [0 -1; -2 3];
lb(:,:,1,3) = [1 1; -1 0];
lb(:,:,2,1) = [-3 2; 0 1];
ub = [];
ub(:,:,1,1) = [1.5 4; 4 10];
ub(:,:,1,2) = [1 2; 0 4];
ub(:,:,1,3) = [2 3; -0.5 2];
ub(:,:,2,1) = [-1 3; 0 2];
I = interval(lb,ub);
Iabs = abs(I);
assert(all(dim(I) == dim(Iabs)));
assert(all(Iabs.inf >= 0,"all"))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
