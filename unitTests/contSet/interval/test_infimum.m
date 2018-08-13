function res = test_infimum(~)
% test_abs - unit_test_function of the infimum of an interval
%
% Syntax:  
%    res = test_infimum( ~ )
%
% Inputs:
%    no
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
% Defenition problem
tol = 1e-9;
res = true;

c = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);

if abs( infimum(c(1)) + 5.0 ) > tol
	res = false;
	disp('test_infimum failed');
	return;
end

if abs( infimum(c(2)) + 4.0 ) > tol
	res = false;
	disp('test_infimum failed');
	return;
end

if abs( infimum(c(3)) + 3.0 ) > tol
	res = false;
	disp('test_infimum failed');
	return;
end

if abs( infimum(c(4)) + 0.0 ) > tol
	res = false;
	disp('test_infimum failed');
	return;
end

if abs( infimum(c(5)) - 0.0 ) > tol
	res = false;
	disp('test_infimum failed');
	return;
end

if abs( infimum(c(6)) - 5.0 ) > tol
	res = false;
	disp('test_infimum failed');
	return;
end

disp('test_infimum successful');
return;

%------------- END OF CODE --------------