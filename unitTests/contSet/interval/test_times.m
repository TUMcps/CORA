function res = test_times( ~ )
% test_times - unit_test_function of times - Overloaded '.*' operator for intervals
%
% Syntax:  
%    res = test_times( ~ )
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
% Written:      04-January-2016
% Last update:  13-January-2016 (DG)
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

a = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
b = interval([-6.1, -4.5, -3.3, 0, 0, 5], [-2.2, 0.0, 2.8, 0, 5.7, 8.2]);
c = a .* b;

if abs( infimum(c(1)) - 4.4 ) > tol || abs( supremum(c(1)) - 30.5 ) > tol
	res = false;
	disp('test_times failed');
	return;
end

if abs( infimum(c(2)) - 0.0 ) > tol || abs( supremum(c(2)) - 18.0 ) > tol
	res = false;
	disp('test_times failed');
	return;
end

if abs( infimum(c(3)) + 8.4 ) > tol || abs( supremum(c(3)) - 9.9 ) > tol
	res = false;
	disp('test_times failed');
	return;
end

if abs( infimum(c(4)) - 0.0 ) > tol || abs( supremum(c(4)) - 0.0 ) > tol
	res = false;
	disp('test_times failed');
	return;
end

if abs( infimum(c(5)) - 0.0 ) > tol || abs( supremum(c(5)) - 28.5 ) > tol
	res = false;
	disp('test_times failed');
	return;
end

if abs( infimum(c(6)) - 25.0 ) > tol || abs( supremum(c(6)) - 65.6 ) > tol
	res = false;
	disp('test_times failed');
	return;
end

disp('test_times successful');
return;

%------------- END OF CODE --------------