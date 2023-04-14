function res = test_interval_isscalar
% test_interval_isscalar - unit test function of isscalar
%
% Syntax:  
%    res = test_interval_isscalar
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

% Author:       Dmitry Grebenyuk
% Written:      16-January-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% define problem
res = true;

test = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = isscalar( test );

if c ~= false
	res = false;
	return;
end

test = interval();
c = isscalar( test );

if c ~= false
	res = false;
	return;
end

test = interval(-5.0, 2);
c = isscalar( test );

if c ~= true
	res = false;
	return;
end

%------------- END OF CODE --------------