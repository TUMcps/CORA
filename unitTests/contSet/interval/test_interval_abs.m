function res = test_interval_abs
% test_interval_abs - unit test function of absolute value
%
% Syntax:  
%    res = test_interval_abs
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
% Last update:  12-January-2016 (DG)
% Last revision:---

%------------- BEGIN CODE --------------

% define problem
tol = 1e-9;
res = true;

test = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = abs( test );

if abs( infimum(c(1)) - 2.0 ) > tol || abs( supremum(c(1)) - 5.0 ) > tol
	res = false;
	disp('test_abs failed');
	return;
end

if abs( infimum(c(2)) - 0.0 ) > tol || abs( supremum(c(2)) - 4.0 ) > tol
	res = false;
	disp('test_abs failed');
	return;
end

if abs( infimum(c(3)) - 0.0 ) > tol || abs( supremum(c(3)) - 3.0 ) > tol
	res = false;
	disp('test_abs failed');
	return;
end

if abs( infimum(c(4)) - 0.0 ) > tol || abs( supremum(c(4)) - 0.0 ) > tol
	res = false;
	disp('test_abs failed');
	return;
end

if abs( infimum(c(5)) - 0.0 ) > tol || abs( supremum(c(5)) - 5.0 ) > tol
	res = false;
	disp('test_abs failed');
	return;
end

if abs( infimum(c(6)) - 5.0 ) > tol || abs( supremum(c(6)) - 8.0 ) > tol
	res = false;
	disp('test_abs failed');
	return;
end

disp('test_abs successful');
return;

%------------- END OF CODE --------------