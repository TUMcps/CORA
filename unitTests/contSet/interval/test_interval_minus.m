function res = test_interval_minus
% test_interval_minus - unit_test_function of minus,
%    overloaded '-' operator for intervals
%
% Syntax:
%    res = test_interval_minus
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

% Authors:       Dmitry Grebenyuk
% Written:       04-January-2016
% Last update:   13-January-2016 (DG)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% scalar tests ---

% equal
a = interval(0, 0);
b = interval(1, 1);
c = a - b;
assert( abs( infimum(c) + 1.0 ) <= tol && abs( supremum(c) + 1.0 ) <= tol )

a = interval(0, 0);
b = interval(1, 1);
c = b - a;
assert( abs( infimum(c) - 1.0 ) <= tol && abs( supremum(c) - 1.0 ) <= tol )

% one negative
a = interval(-1, 0);
b = interval(-1, 0);
c = a - b;
assert( abs( infimum(c) + 1.0 ) <= tol && abs( supremum(c) - 1.0 ) <= tol )

a = interval(-1, 0);
b = interval(-2, 0);
c = a - b;
assert( abs( infimum(c) + 1.0 ) <= tol && abs( supremum(c) - 2.0 ) <= tol )

a = interval(-1.0, 0);
b = interval(-2.0, 0);
c = a - b;
assert( abs( infimum(c) + 1.0 ) <= tol && abs( supremum(c) - 2.0 ) <= tol )

% both positive
a = interval(0, 1);
b = interval(0, 1);
c = a - b;
assert( abs( infimum(c) + 1.0 ) <= tol && abs( supremum(c) - 1.0 ) <= tol )

a = interval(0, 1);
b = interval(0, 2);
c = a - b;
assert( abs( infimum(c) + 2.0 ) <= tol && abs( supremum(c) - 1.0 ) <= tol )

a = interval(0, 1.0);
b = interval(0, 2.0);
c = a - b;
assert( abs( infimum(c) + 2.0 ) <= tol && abs( supremum(c) - 1.0 ) <= tol )

% mixed
a = interval(-2.0, 2.0);
b = interval(-3.0, 2.0);
c = a - b;
assert( abs( infimum(c) + 4.0 ) <= tol && abs( supremum(c) - 5.0 ) <= tol )

a = interval(-2.0, 2.0);
b = interval(-3.0, 2.0);
c = b - a;
assert( abs( infimum(c) + 5.0 ) <= tol && abs( supremum(c) - 4.0 ) <= tol )

a = interval(-2.0, 1.0);
c = 1 - a;
assert( abs( infimum(c) + 0.0 ) <= tol && abs( supremum(c) - 3.0 ) <= tol )

a = interval(-2.0, 1.0);
c = a - 1;
assert( abs( infimum(c) + 3.0 ) <= tol && abs( supremum(c) - 0.0 ) <= tol )

% n-d test
a = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
b = interval([-6.1, -4.5, -3.3, 0, 0, 5], [-2.2, 0.0, 2.8, 0, 5.7, 8.2]);
c = a - b;
assert( abs( infimum(c(1)) + 2.8 ) <= tol && abs( supremum(c(1)) - 4.1 ) <= tol )
assert( abs( infimum(c(2)) + 4.0 ) <= tol && abs( supremum(c(2)) - 4.5 ) <= tol )
assert( abs( infimum(c(3)) + 5.8 ) <= tol && abs( supremum(c(3)) - 5.3 ) <= tol )
assert( abs( infimum(c(4)) + 0.0 ) <= tol && abs( supremum(c(4)) - 0.0 ) <= tol )
assert( abs( infimum(c(5)) + 5.7 ) <= tol && abs( supremum(c(5)) - 5.0 ) <= tol )
assert( abs( infimum(c(6)) + 3.2 ) <= tol && abs( supremum(c(6)) - 3.0 ) <= tol )

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
