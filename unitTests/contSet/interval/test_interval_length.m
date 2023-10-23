function res = test_interval_length
% test_interval_length - unit_test_function of length
%
% Syntax:
%    res = test_interval_length
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

% Authors:       Dmitry Grebenyuk
% Written:       19-January-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true;

a = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = length(a);

if c ~= 6
	res = false;
	return;
end

a = interval([-5.0; -4.0; -3; 0; 0; 5], [-2; 0.0; 2.0; 0; 5; 8]);
c = length(a);

if c ~= 6
	res = false;
	return;
end

a = interval();
c = length(a);

if c ~= 0
	res = false;
	return;
end

% ------------------------------ END OF CODE ------------------------------
