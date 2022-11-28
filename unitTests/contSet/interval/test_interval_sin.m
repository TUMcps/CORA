function res = test_interval_sin
% test_interval_sin - unit_test_function of sine for intervals,
%    overloaded 'sin()' function for intervals
%
% Syntax:  
%    res = test_interval_sin
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

% Author:       Dmitry Grebenyuk
% Written:      13-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;


a = interval([-0; 0.0; 0.0; 0; 0; -pi/4; -2*pi], [pi/4; pi/2; pi; 2*pi; 4*pi; pi/4; 4*pi]);
c = sin(a);

if abs( infimum(c(1)) - 0.0 ) > tol || abs( supremum(c(1)) - 0.707106781 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2)) - 0.0 ) > tol || abs( supremum(c(2)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3)) - 0.0 ) > tol || abs( supremum(c(3)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(4)) + 1.0 ) > tol || abs( supremum(c(4)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(5)) + 1.0 ) > tol || abs( supremum(c(5)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(6)) + 0.707106781 ) > tol || abs( supremum(c(6)) - 0.707106781 ) > tol
	res = false;
	return;
end

if abs( infimum(c(7)) + 1.0 ) > tol || abs( supremum(c(7)) - 1.0 ) > tol
	res = false;
	return;
end


a = interval([-0, 0.0, 0.0, 0, 0, -pi/4, -2*pi], [pi/4, pi/2, pi, 2*pi, 4*pi, pi/4, 4*pi]);
c = sin(a);

if abs( infimum(c(1)) - 0.0 ) > tol || abs( supremum(c(1)) - 0.707106781 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2)) - 0.0 ) > tol || abs( supremum(c(2)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3)) - 0.0 ) > tol || abs( supremum(c(3)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(4)) + 1.0 ) > tol || abs( supremum(c(4)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(5)) + 1.0 ) > tol || abs( supremum(c(5)) - 1.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(6)) + 0.707106781 ) > tol || abs( supremum(c(6)) - 0.707106781 ) > tol
	res = false;
	return;
end

if abs( infimum(c(7)) + 1.0 ) > tol || abs( supremum(c(7)) - 1.0 ) > tol
	res = false;
	return;
end

%------------- END OF CODE --------------