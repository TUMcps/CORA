function res = test_and( ~ )

% test_and - unit_test_function of logical conjunction - Overloaded '&' operator for intervals
%
% Syntax:  
%    res = test_and( ~ )
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
% Written:      05-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

a = interval([-10; -2.0; 10.0], [-5; 8.0; 15.0]);
b = interval([-11; 2; 11.0], [-6; 9; 12.0]);
c = a & b;
if abs( infimum(c(1)) + 10.0 ) > tol || abs( supremum(c(1)) + 6.0 ) > tol
	res = false;
	disp('test_and failed');
	return;
end

if abs( infimum(c(2)) - 2.0 ) > tol || abs( supremum(c(2)) - 8.0 ) > tol
	res = false;
	disp('test_and failed');
	return;
end

if abs( infimum(c(3)) - 11.0 ) > tol || abs( supremum(c(3)) - 12.0 ) > tol
	res = false;
	disp('test_and failed');
	return;
end

disp('test_and successful');
return;

%------------- END OF CODE --------------