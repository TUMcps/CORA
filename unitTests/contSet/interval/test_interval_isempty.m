function res = test_interval_isempty
% test_interval_isempty - unit_test_function of isempty
%
% Syntax:  
%    res = test_interval_isempty
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

% Author:       Dmitry Grebenyuk
% Written:      16-January-2016
% Last update:  23-March-2021 (MW, rewrite syntax)
% Last revision:---

%------------- BEGIN CODE --------------

% define problem
res = true;

Int = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
if isempty(Int)
	res = false;
	disp('test_isempty failed');
	return;
end

Int = interval();
if ~isempty(Int)
	res = false;
	disp('test_isempty failed');
	return;
end

Int = interval(-5.0, 2);
if isempty(Int)
	res = false;
	disp('test_isempty failed');
	return;
end

disp('test_isempty successful');

%------------- END OF CODE --------------