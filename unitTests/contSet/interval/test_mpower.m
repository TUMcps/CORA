function res = test_mpower( ~ )

% test_mpower - unit_test_function of power - Overloaded '^' operator for intervals
%
% Syntax:  
%    res = test_mpower( ~ )
%
% Inputs:
%    intVal - no
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

% Author:       Dmitry Grebenyuk
% Written:      05-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

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
