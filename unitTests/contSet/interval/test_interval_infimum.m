function res = test_interval_infimum
% test_interval_infimum - unit test function of the infimum of an interval
%
% Syntax:  
%    res = test_interval_infimum
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

tol = 1e-9;
res = true;

c = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);

if abs( infimum(c(1)) + 5.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(2)) + 4.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(3)) + 3.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(4)) + 0.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(5)) - 0.0 ) > tol
	res = false;
	return;
end

if abs( infimum(c(6)) - 5.0 ) > tol
	res = false;
	return;
end

%------------- END OF CODE --------------