function res = test_transpose( ~ )

% test_transpose - unit_test_function of transpose - Overloaded '.'' operator for intervals
%
% Syntax:  
%    res = test_transpose( ~ )
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
%
% Author:       Dmitry Grebenyuk
% Written:      07-February-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

a = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
b = a + 1;
c = b + 2;
d = [a; b; c];
c = d.';


if abs( infimum(c(1, 1)) + 5.0 ) > tol || abs( supremum(c(1, 1)) + 2.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(2, 1)) + 4.0 ) > tol || abs( supremum(c(2, 1)) - 0.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(3, 1)) + 3.0 ) > tol || abs( supremum(c(3, 1)) - 2.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(4, 1)) - 0.0 ) > tol || abs( supremum(c(4, 1)) + 0.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(5, 1)) + 0.0 ) > tol || abs( supremum(c(5, 1)) - 5.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(6, 1)) - 5.0 ) > tol || abs( supremum(c(6, 1)) - 8.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(1, 2)) + 4.0 ) > tol || abs( supremum(c(1, 2)) + 1.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(2, 2)) + 3.0 ) > tol || abs( supremum(c(2, 2)) - 1.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(3, 2)) + 2.0 ) > tol || abs( supremum(c(3, 2)) - 3.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(4, 2)) - 1.0 ) > tol || abs( supremum(c(4, 2)) - 1.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(5, 2)) - 1.0 ) > tol || abs( supremum(c(5, 2)) - 6.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(6, 2)) - 6.0 ) > tol || abs( supremum(c(6, 2)) - 9.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(1, 3)) + 2.0 ) > tol || abs( supremum(c(1, 3)) - 1.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(2, 3)) + 1.0 ) > tol || abs( supremum(c(2, 3)) - 3.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(3, 3)) + 0.0 ) > tol || abs( supremum(c(3, 3)) - 5.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(4, 3)) - 3.0 ) > tol || abs( supremum(c(4, 3)) - 3.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(5, 3)) - 3.0 ) > tol || abs( supremum(c(5, 3)) - 8.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

if abs( infimum(c(6, 3)) - 8.0 ) > tol || abs( supremum(c(6, 3)) - 11.0 ) > tol
	res = false;
	disp('test_transpose failed');
	return;
end

disp('test_transpose successful');
return;

%------------- END OF CODE --------------
