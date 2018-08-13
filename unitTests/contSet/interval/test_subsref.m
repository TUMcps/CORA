function res = test_subsref( ~ )

% test_subsref - unit_test_function of subscripted reference for objects
%
% Syntax:  
%    res = test_subsref( ~ )
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
% Written:      19-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

a = interval([-5.0; -4.0; -3; 0; 0; 5], [-2; 0.0; 2.0; 0; 5; 8]);
c = a;

if abs( infimum(c(1,1)) + 5.0 ) > tol || abs( supremum(c(1,1)) + 2.0 ) > tol
	res = false;
	disp('test_subsref failed');
	return;
end

if abs( infimum(c(2,1)) + 4.0 ) > tol || abs( supremum(c(2,1)) + 0.0 ) > tol
	res = false;
	disp('test_subsref failed');
	return;
end

if abs( infimum(c(3,1)) + 3.0 ) > tol || abs( supremum(c(3,1)) - 2.0 ) > tol
	res = false;
	disp('test_subsref failed');
	return;
end

if abs( infimum(c(4,1)) - 0.0 ) > tol || abs( supremum(c(4,1)) - 0.0 ) > tol
	res = false;
	disp('test_subsref failed');
	return;
end

if abs( infimum(c(5,1)) - 0.0 ) > tol || abs( supremum(c(5,1)) - 5.0 ) > tol
	res = false;
	disp('test_subsref failed');
	return;
end

if abs( infimum(c(6,1)) - 5.0 ) > tol || abs( supremum(c(6,1)) - 8.0 ) > tol
	res = false;
	disp('test_subsref failed');
	return;
end

disp('test_subsref successful');
return;

%------------- END OF CODE --------------