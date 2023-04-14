function res = test_interval_horzcat
% test_interval_horzcat - unit test function of the operator for
%    horizontal concatenation, e.g. a = [b,c,d];
%
% Syntax:  
%    res = test_interval_horzcat
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

% Author:       Dmitry Grebenyuk
% Written:      16-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% define problem
tol = 1e-9;
res = true;

a = interval([-5.0; -4.0; -3; 0; 0; 5], [-2; 0.0; 2.0; 0; 5; 8]);
b = a + 1;
c = b + 2;
c = [a, b, c];

if abs( infimum(c(1, 1)) + 5.0 ) > tol || abs( supremum(c(1, 1)) + 2.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 1)) + 4.0 ) > tol || abs( supremum(c(2, 1)) - 0.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 1)) + 3.0 ) > tol || abs( supremum(c(3, 1)) - 2.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(4, 1)) - 0.0 ) > tol || abs( supremum(c(4, 1)) + 0.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(5, 1)) + 0.0 ) > tol || abs( supremum(c(5, 1)) - 5.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(6, 1)) - 5.0 ) > tol || abs( supremum(c(6, 1)) - 8.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(1, 2)) + 4.0 ) > tol || abs( supremum(c(1, 2)) + 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 2)) + 3.0 ) > tol || abs( supremum(c(2, 2)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 2)) + 2.0 ) > tol || abs( supremum(c(3, 2)) - 3.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(4, 2)) - 1.0 ) > tol || abs( supremum(c(4, 2)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(5, 2)) - 1.0 ) > tol || abs( supremum(c(5, 2)) - 6.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(6, 2)) - 6.0 ) > tol || abs( supremum(c(6, 2)) - 9.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(1, 3)) + 2.0 ) > tol || abs( supremum(c(1, 3)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2, 3)) + 1.0 ) > tol || abs( supremum(c(2, 3)) - 3.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3, 3)) + 0.0 ) > tol || abs( supremum(c(3, 3)) - 5.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(4, 3)) - 3.0 ) > tol || abs( supremum(c(4, 3)) - 3.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(5, 3)) - 3.0 ) > tol || abs( supremum(c(5, 3)) - 8.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(6, 3)) - 8.0 ) > tol || abs( supremum(c(6, 3)) - 11.0 ) > tol
	res = false;
	return;
end

%------------- END OF CODE --------------