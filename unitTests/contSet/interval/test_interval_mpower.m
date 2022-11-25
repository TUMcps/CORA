function res = test_interval_mpower
% test_interval_mpower - unit_test_function of power,
%    overloaded '^' operator for intervals
%
% Syntax:  
%    res = test_interval_mpower
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:      05-January-2016
% Last update:  08-August-2020 (MW, extend by random tests)
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

% 1. analytical tests -----------------------------------------------------

a = interval(0, 2);
c = a ^ 1;
if abs( infimum(c) - 0.0 ) > tol || abs( supremum(c) - 2.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(0, 2);
c = a ^ 2;
if abs( infimum(c) - 0.0 ) > tol || abs( supremum(c) - 4.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(-2, 0);
c = a ^ 2;
if abs( infimum(c) - 0.0 ) > tol || abs( supremum(c) - 4.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(-2, 0);
c = a ^ 3;
if abs( infimum(c) + 8.0 ) > tol || abs( supremum(c) - 0.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(-3, 2);
c = a ^ 2;
if abs( infimum(c) - 0.0 ) > tol || abs( supremum(c) - 9.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(-3, 2);
c = a ^ 3;
if abs( infimum(c) + 27.0 ) > tol || abs( supremum(c) - 8.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(-3, -2);
c = a ^ 2;
if abs( infimum(c) - 4.0 ) > tol || abs( supremum(c) - 9.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(-3, -2);
c = a ^ 3;
if abs( infimum(c) + 27.0 ) > tol || abs( supremum(c) + 8.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(2, 3);
c = a ^ 2;
if abs( infimum(c) - 4.0 ) > tol || abs( supremum(c) - 9.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end

a = interval(2, 3);
c = a ^ 3;
if abs( infimum(c) - 8.0 ) > tol || abs( supremum(c) - 27.0 ) > tol
	res = false;
	disp('test_mpower failed');
	return;
end


disp('test_mpower successful');
return;

%------------- END OF CODE --------------
