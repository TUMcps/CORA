function res = test_interval_asin
% test_interval_asin - unit test function of sine for intervals,
%    overloaded 'asin()' function for intervals
%
% Syntax:
%    res = test_interval_asin
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

% Authors:       Mark Wetzlinger
% Written:       08-August-2020
% Last update:   03-December-2023 (MW, add out of bounds cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-9;

% Check special values
% x     asin(x)
% -1    -pi/2
% 0     0
% 1     pi/2

I = interval([-1,0],[1,0]);
I_asin = asin(I);
I_true = interval([-pi/2,0],[pi/2,0]);
res(end+1,1) = isequal(I_asin,I_true,tol);

% check out of bounds
I = interval(-2,0);
try
    I_asin = acos(I);
    res(end+1,1) = false;
end
I = interval(0.5,1.1);
try
    I_asin = acos(I);
    res(end+1,1) = false;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
