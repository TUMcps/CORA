function res = test_interval_vertcat
% test_interval_vertcat - unit test function of the operator for
%    vertical concatenation, e.g. a = [b;c;d];
%
% Syntax:
%    res = test_interval_vertcat
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
% Written:       19-January-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define problem
tol = 1e-9;
res = true;

a = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
b = a + 1;
c = b + 2;
c = [a; b; c];

if abs( infimum(c(1, 1)) + 5.0 ) > tol || abs( supremum(c(1, 1)) + 2.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(1, 2)) + 4.0 ) > tol || abs( supremum(c(1, 2)) - 0.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(1, 3)) + 3.0 ) > tol || abs( supremum(c(1, 3)) - 2.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(1, 4)) - 0.0 ) > tol || abs( supremum(c(1, 4)) + 0.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(1, 5)) + 0.0 ) > tol || abs( supremum(c(1, 5)) - 5.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(1, 6)) - 5.0 ) > tol || abs( supremum(c(1, 6)) - 8.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 1)) + 4.0 ) > tol || abs( supremum(c(2, 1)) + 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 2)) + 3.0 ) > tol || abs( supremum(c(2, 2)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 3)) + 2.0 ) > tol || abs( supremum(c(2, 3)) - 3.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 4)) - 1.0 ) > tol || abs( supremum(c(2, 4)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 5)) - 1.0 ) > tol || abs( supremum(c(2, 5)) - 6.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 6)) - 6.0 ) > tol || abs( supremum(c(2, 6)) - 9.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 1)) + 2.0 ) > tol || abs( supremum(c(3, 1)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 2)) + 1.0 ) > tol || abs( supremum(c(3, 2)) - 3.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 3)) + 0.0 ) > tol || abs( supremum(c(3, 3)) - 5.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 4)) - 3.0 ) > tol || abs( supremum(c(3, 4)) - 3.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 5)) - 3.0 ) > tol || abs( supremum(c(3, 5)) - 8.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 6)) - 8.0 ) > tol || abs( supremum(c(3, 6)) - 11.0 ) > tol
	res = false;
	return;
end

% ------------------------------ END OF CODE ------------------------------
