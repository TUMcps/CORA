function res = test_interval_volume
% test_interval_volume - unit test function of volume
%
% Syntax:
%    res = test_interval_volume
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
% Written:       28-August-2019
% Last update:   04-December-2023 (MW, add degenerate and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-9;

% empty
I = interval.empty(2);
vol = volume(I);
vol_true = 0;
res(end+1,1) = withinTol(vol,vol_true,tol);

% bounded, full-dimensional
I = interval([-2; -4; -3],[3; 1; 2]);
vol = volume(I);
vol_true = 125;
res(end+1,1) = withinTol(vol,vol_true,tol);

% bounded, degenerate
I = interval([0;2],[1;2]);
vol = volume(I);
vol_true = 0;
res(end+1,1) = withinTol(vol,vol_true,tol);

% unbounded, full-dimensional
I = interval([-Inf;-2],[1;2]);
vol = volume(I);
vol_true = Inf;
res(end+1,1) = withinTol(vol,vol_true,tol);

% unbounded, degenerate
I = interval([-Inf;2],[1;2]);
vol = volume(I);
vol_true = 0;
res(end+1,1) = withinTol(vol,vol_true,tol);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
