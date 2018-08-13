function res = test_log( ~ )

% test_log - unit_test_function of natural logarithm for intervals - Overloaded 'log()' function for intervals
%
% Syntax:  
%    res = test_log( ~ )
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
%
% Author:       Dmitry Grebenyuk
% Written:      07-February-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

a = interval([-2, -2; 0, 1], [-1, 0; 2, 2]);
c = log(a);

if isnan(infimum(c(1, 1))) ~= 1 || isnan(supremum(c(1, 1))) ~= 1
	res = false;
	disp('test_log failed');
	return;
end

if isnan(infimum(c(1, 2))) ~= 1 || isinf(supremum(c(1, 2))) ~= 1
	res = false;
	disp('test_log failed');
	return;
end

if isinf(infimum(c(2, 1))) ~= 1 || abs( supremum(c(2, 1)) - 0.6931471805599 ) > tol
	res = false;
	disp('test_log failed');
	return;
end

if abs( infimum(c(2, 2)) - 0.0 ) > tol || abs( supremum(c(2, 2)) - 0.6931471805599 ) > tol
	res = false;
	disp('test_log failed');
	return;
end



disp('test_log successful');
return;

%------------- END OF CODE --------------