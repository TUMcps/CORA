function res = test_interval_supremum
% test_interval_supremum - unit test function of the supremum of an interval
%
% Syntax:  
%    res = test_interval_supremum
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
% Written:      14-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% define problem
tol = 1e-9;
res = true;

c = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);

if abs( supremum(c(1)) + 2.0 ) > tol
	res = false;
	disp('test_supremum failed');
	return;
end

if abs( supremum(c(2)) + 0.0 ) > tol
	res = false;
	disp('test_supremum failed');
	return;
end

if abs( supremum(c(3)) - 2.0 ) > tol
	res = false;
	disp('test_supremum failed');
	return;
end

if abs( supremum(c(4)) - 0.0 ) > tol
	res = false;
	disp('test_supremum failed');
	return;
end

if abs( supremum(c(5)) - 5.0 ) > tol
	res = false;
	disp('test_supremum failed');
	return;
end

if abs( supremum(c(6)) - 8.0 ) > tol
	res = false;
	disp('test_supremum failed');
	return;
end

disp('test_supremum successful');
return;

%------------- END OF CODE --------------
