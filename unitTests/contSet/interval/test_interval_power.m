function res = test_interval_power
% test_interval_power - unit_test_function of power,
%    overloaded '.^' operator for intervals
%
% Syntax:
%    res = test_interval_power
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

% Authors:       Tobias Ladner
% Written:       14-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% bounded
I = interval(0, 2);
I_power = I .^ 1;
assert(isequal(I_power,I,tol));

I = interval(0, 2);
I_power = I .^ 2;
I_true = interval(0,4);
assert(isequal(I_power,I_true,tol));

% negative
I = interval(-2, 0);
I_power = I .^ 2;
I_true = interval(0,4);
assert(isequal(I_power,I_true,tol));

I = interval(-2, 0);
I_power = I .^ 3;
I_true = interval(-8,0);
assert(isequal(I_power,I_true,tol));

I = interval(-3, 2);
I_power = I ^ 2;
I_true = interval(0,9);
assert(isequal(I_power,I_true,tol));

I = interval(-3, 2);
I_power = I ^ 3;
I_true = interval(-27,8);
assert(isequal(I_power,I_true,tol));

% both negative
I = interval(-3, -2);
I_power = I .^ 2;
I_true = interval(4,9);
assert(isequal(I_power,I_true,tol));

I = interval(-3, -2);
I_power = I .^ 3;
I_true = interval(-27,-8);
assert(isequal(I_power,I_true,tol));

% both positive
I = interval(2, 3);
I_power = I .^ 2;
I_true = interval(4,9);
assert(isequal(I_power,I_true,tol));

I = interval(2, 3);
I_power = I .^ 3;
I_true = interval(8,27);
assert(isequal(I_power,I_true,tol));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
