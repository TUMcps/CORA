function res = test_length( ~ )

% test_length - unit_test_function of length()
%
% Syntax:  
%    res = test_length( ~ )
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
c = length(a);

if c ~= 6
	res = false;
	disp('test_length failed');
	return;
end

a = interval([-5.0; -4.0; -3; 0; 0; 5], [-2; 0.0; 2.0; 0; 5; 8]);
c = length(a);

if c ~= 6
	res = false;
	disp('test_length failed');
	return;
end

a = interval();
c = length(a);

if c ~= 0
	res = false;
	disp('test_length failed');
	return;
end

disp('test_length successful');
return;

%------------- END OF CODE --------------