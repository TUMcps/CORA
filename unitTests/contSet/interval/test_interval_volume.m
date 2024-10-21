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

tol = 1e-9;

% empty
I = interval.empty(2);
vol = volume(I);
vol_true = 0;
assert(withinTol(vol,vol_true,tol));

% bounded, full-dimensional
I = interval([-2; -4; -3],[3; 1; 2]);
vol = volume(I);
vol_true = 125;
assert(withinTol(vol,vol_true,tol));

% bounded, degenerate
I = interval([0;2],[1;2]);
vol = volume(I);
vol_true = 0;
assert(withinTol(vol,vol_true,tol));

% unbounded, full-dimensional
I = interval([-Inf;-2],[1;2]);
vol = volume(I);
vol_true = Inf;
assert(withinTol(vol,vol_true,tol));

% unbounded, degenerate
I = interval([-Inf;2],[1;2]);
vol = volume(I);
vol_true = 0;
assert(withinTol(vol,vol_true,tol));

% matrix
lb = [ -0.265 -0.520 ; -0.713 -1.349 ];
ub = [ 2.496 0.546 ; 1.357 4.379 ];
I = interval(lb,ub);
vol = volume(I);
vol_true = 34.8977;
assert(withinTol(vol,vol_true,1e-3));

% n-d arrays
lb = reshape([ -0.785 -2.226 -0.776 -1.997 -3.416 -2.200 0.121 -0.638 -0.458 -3.367 -2.707 -3.100 ], [2,2,3]);
ub = reshape([ 2.633 -1.267 0.718 -0.164 2.115 2.613 1.828 -0.064 0.519 0.399 0.426 1.443 ], [2,2,3]);
I = interval(lb,ub);
vol = volume(I);
vol_true = 12261.587;
assert(withinTol(vol,vol_true,1e-2));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
