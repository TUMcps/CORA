function res = test_exp(~)
% test_exp - unit_test_function of exponent
%
% Syntax:  
%    res = test_exp( ~ )
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

test = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = exp( test );

if abs( infimum(c(1)) - 0.006737947 ) > tol || abs( supremum(c(1)) - 0.135335283 ) > tol
	res = false;
	disp('test_exp failed');
	return;
end

if abs( infimum(c(2)) - 0.0183156389 ) > tol || abs( supremum(c(2)) - 1.0 ) > tol
	res = false;
	disp('test_exp failed');
	return;
end

if abs( infimum(c(3)) - 0.0497870684 ) > tol || abs( supremum(c(3)) - 7.3890560989 ) > tol
	res = false;
	disp('test_exp failed');
	return;
end

if abs( infimum(c(4)) - 1.0 ) > tol || abs( supremum(c(4)) - 1.0 ) > tol
	res = false;
	disp('test_exp failed');
	return;
end

if abs( infimum(c(5)) - 1.0 ) > tol || abs( supremum(c(5)) - 148.4131591026 ) > tol
	res = false;
	disp('test_exp failed');
	return;
end

if abs( infimum(c(6)) - 148.4131591026 ) > tol || abs( supremum(c(6)) - 2980.9579870418 ) > tol
	res = false;
	disp('test_exp failed');
	return;
end

disp('test_exp successful');
return;

%------------- END OF CODE --------------