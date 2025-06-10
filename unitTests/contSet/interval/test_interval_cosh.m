function res = test_interval_cosh
% test_interval_cosh - unit test function of cosh for intervals,
%    overloaded 'cosh()' function for intervals
%
% Syntax:
%    res = test_interval_cosh
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
I = interval(-1,2);
I_cosh = cosh(I);
I_true = interval(1,3.762195691083631);
assert(isequal(I_cosh,I_true,tol));

% unbounded
I = interval(-inf,2);
I_cosh = cosh(I);
I_true = interval(1,Inf);
assert(isequal(I_cosh,I_true,tol));

I = interval(-1,inf);
I_cosh = cosh(I);
I_true = interval(1,Inf);
assert(isequal(I_cosh,I_true,tol));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
