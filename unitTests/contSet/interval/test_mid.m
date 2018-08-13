function res = test_mid( ~ )

% test_mid - unit_test_function of mid()
%
% Syntax:  
%    res = test_mid( ~ )
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
% Written:      15-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

a = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = mid(a);

if abs( c(1) + 3.5 ) > tol
	res = false;
	disp('test_mid failed');
	return;
end

if abs( c(2) + 2 ) > tol
	res = false;
	disp('test_mid failed');
	return;
end

if abs( c(3) + 0.5 ) > tol
	res = false;
	disp('test_mid failed');
	return;
end

if abs( c(4) + 0 ) > tol
	res = false;
	disp('test_mid failed');
	return;
end

if abs( c(5) - 2.5 ) > tol
	res = false;
	disp('test_mid failed');
	return;
end

if abs( c(6) - 6.5 ) > tol
	res = false;
	disp('test_mid failed');
	return;
end

%%--------

a = interval();
c = mid(a);

if isempty(c) ~= true
	res = false;
	disp('test_mid failed');
	return;
end

%%--------

disp('test_mid successful');
return;

%------------- END OF CODE --------------