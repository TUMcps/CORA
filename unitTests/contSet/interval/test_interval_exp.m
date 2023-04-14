function res = test_interval_exp
% test_interval_exp - unit test function of exponential function
%
% Syntax:  
%    res = test_interval_exp
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
% Written:      14-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% define problem
tol = 1e-9;
res = true;

test = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = exp( test );

if abs( infimum(c(1)) - 0.006737947 ) > tol || abs( supremum(c(1)) - 0.135335283 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2)) - 0.0183156389 ) > tol || abs( supremum(c(2)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3)) - 0.0497870684 ) > tol || abs( supremum(c(3)) - 7.3890560989 ) > tol
	res = false;
	return;
end

if abs( infimum(c(4)) - 1.0 ) > tol || abs( supremum(c(4)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(5)) - 1.0 ) > tol || abs( supremum(c(5)) - 148.4131591026 ) > tol
	res = false;
	return;
end

if abs( infimum(c(6)) - 148.4131591026 ) > tol || abs( supremum(c(6)) - 2980.9579870418 ) > tol
	res = false;
	return;
end

%------------- END OF CODE --------------