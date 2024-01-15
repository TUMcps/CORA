function res = test_interval_sqrt
% test_interval_sqrt - unit test function of sqrt
%
% Syntax:
%    res = test_interval_sqrt
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
% Written:       29-August-2019
% Last update:   04-December-2023 (MW, add empty and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-9;

% empty
I = interval.empty(2);
I_sqrt = sqrt(I);
res(end+1,1) = representsa(I_sqrt,'emptySet');

% bounded (only perfect squares)
I = interval([4; 9; 4; 16; 1], [9; 25; 36; 100; 4]);
I_sqrt = sqrt(I);
I_true = interval([2; 3; 2; 4; 1], [3; 5; 6; 10; 2]);
res(end+1,1) = isequal(I_sqrt,I_true,tol);

% unbounded
I = interval([2;4],[Inf;9]);
I_sqrt = sqrt(I);
I_true = interval([sqrt(2);2],[Inf;3]);
res(end+1,1) = isequal(I_sqrt,I_true,tol);

% out of bounds
I = interval(-2,1);
try
    sqrt(I);
    res = false;
end


% add results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
