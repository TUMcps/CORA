function res = test_isempty(~)
% test_isempty - unit_test_function of isempty
%
% Syntax:  
%    res = test_isempty( ~ )
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% Defenition problem
res = true;

test = interval([-5.0, -4.0, -3, 0, 0, 5], [-2, 0.0, 2.0, 0, 5, 8]);
c = isempty( test );

if c ~= false
	res = false;
	disp('test_isempty failed');
	return;
end

test = interval();
c = isempty( test );

if c ~= true
	res = false;
	disp('test_isempty failed');
	return;
end

test = interval(-5.0, 2);
c = isempty( test );

if c ~= false
	res = false;
	disp('test_isempty failed');
	return;
end

disp('test_isempty successful');
return;

%------------- END OF CODE --------------