function res = test_interval_center
% test_interval_center - unit test function of center
%
% Syntax:  
%    res = test_interval_center
%
% Inputs:
%    .
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author:       Dmitry Grebenyuk
% Written:      15-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

% empty set
if ~isnumeric(center(interval())) || ~isempty(center(interval()))
    res = false;
    return
end

a = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = center(a);

if abs( c(1) + 3.5 ) > tol
	res = false;
	return;
end

if abs( c(2) + 2 ) > tol
	res = false;
	return;
end

if abs( c(3) + 0.5 ) > tol
	res = false;
	return;
end

if abs( c(4) + 0 ) > tol
	res = false;
	return;
end

if abs( c(5) - 2.5 ) > tol
	res = false;
	return;
end

if abs( c(6) - 6.5 ) > tol
	res = false;
	return;
end

%------------- END OF CODE --------------