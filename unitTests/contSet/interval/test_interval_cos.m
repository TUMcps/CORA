function res = test_interval_cos
% test_interval_cos - unit test function of cosine for intervals,
%    overloaded 'cos()' function for intervals
%
% Syntax:  
%    res = test_interval_cos
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
% Written:      05-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;


a = interval([-0; 0.0; 0.0; 0; 0; -pi/4; -2*pi], [pi/4; pi/2; pi; 2*pi; 4*pi; pi/4; 4*pi]);
c = cos(a);

if abs( infimum(c(1)) - 0.707106781 ) > tol || abs( supremum(c(1)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(2)) - 0.0 ) > tol || abs( supremum(c(2)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(3)) + 1.0 ) > tol || abs( supremum(c(3)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(4)) + 1.0 ) > tol || abs( supremum(c(4)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(5)) + 1.0 ) > tol || abs( supremum(c(5)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(6)) - 0.707106781 ) > tol || abs( supremum(c(6)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(7)) + 1.0 ) > tol || abs( supremum(c(7)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

a = interval([-0, 0.0, 0.0, 0, 0, -pi/4, -2*pi], [pi/4, pi/2, pi, 2*pi, 4*pi, pi/4, 4*pi]);
c = cos(a);

if abs( infimum(c(1)) - 0.707106781 ) > tol || abs( supremum(c(1)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(2)) - 0.0 ) > tol || abs( supremum(c(2)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(3)) + 1.0 ) > tol || abs( supremum(c(3)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(4)) + 1.0 ) > tol || abs( supremum(c(4)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(5)) + 1.0 ) > tol || abs( supremum(c(5)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(6)) - 0.707106781 ) > tol || abs( supremum(c(6)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end

if abs( infimum(c(7)) + 1.0 ) > tol || abs( supremum(c(7)) - 1.0 ) > tol
	res = false;
	disp('test_cos failed');
	return;
end


disp('test_cos successful');
return;

%------------- END OF CODE --------------