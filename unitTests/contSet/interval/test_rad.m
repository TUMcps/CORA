function res = test_rad( ~ )

% test_mid - unit_test_function of rad()
%
% Syntax:  
%    res = test_rad( ~ )
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

a = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = rad(a);

if abs( c(1) - 1.5 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(2) - 2.0 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(3) - 2.5 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(4) - 0 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(5) - 2.5 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(6) - 1.5 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

%%-------

a = interval([-5.0; -4.0; -3; 0; 0; 5], [-2; 0.0; 2.0; 0; 5; 8]);
c = rad(a);

if abs( c(1) - 1.5 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(2) - 2.0 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(3) - 2.5 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(4) - 0 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(5) - 2.5 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

if abs( c(6) - 1.5 ) > tol
	res = false;
	disp('test_rad failed');
	return;
end

%%--------

a = interval();
c = rad(a);

if isempty(c) ~= true
	res = false;
	disp('test_rad failed');
	return;
end

%%--------

disp('test_mid successful');
return;

%------------- END OF CODE --------------