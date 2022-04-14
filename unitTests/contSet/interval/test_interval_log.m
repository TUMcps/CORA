function res = test_interval_log
% test_interval_log - unit test function of natural logarithm for intervals
%    overloaded 'log()' function for intervals
%
% Syntax:  
%    res = test_interval_log
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:      07-February-2016
% Last update:  08-June-2020 (MW, rewrite based on new NaN/Inf handling)
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

Int1 = interval([-2; -2],[-1; 0]);

try
    c = log(Int1);
    res = false;
    disp('test_log failed');
    return;
catch
    % log results in NaN -> should throw error
end

Int2 = interval([0, 1],[2, 2]);
c = log(Int2);

if ~isinf(infimum(c(1, 1))) || abs( supremum(c(1, 1)) - 0.6931471805599 ) > tol
	res = false;
	disp('test_log failed');
	return;
end

if abs( infimum(c(1, 2)) - 0.0 ) > tol || abs( supremum(c(1, 2)) - 0.6931471805599 ) > tol
	res = false;
	disp('test_log failed');
	return;
end


disp('test_log successful');

%------------- END OF CODE --------------