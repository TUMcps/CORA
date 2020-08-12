function res = test_containsPoint
% test_containsPoint - unit test function of containsPoint
%
% Syntax:  
%    res = test_containsPoint
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
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create interval
lower = [-3; -9; -4; -7; -1];
upper = [4;   2;  6;  3;  8];
Int = interval(lower, upper);

% points inside
p_inside = [0, 0, 0, 0, 0;
            1,-4, 3,-6, 5;
           -2,-6,-2, 2, 7;
           -3, 2, 6,-7, 8]';
res_inside = all(containsPoint(Int,p_inside));

% points outside
p_outside = [5, 3, 7, 4, 9;
             1,-4,-5,-6, 5;
            -2, 3,-2, 6, 7;
            -3, 2, 6,-7,10]';
res_outside = ~any(containsPoint(Int,p_outside));

% check if all points are in interval
res = res_inside && res_outside;

if res
    disp('test_containsPoint successful');
else
    disp('test_containsPoint failed');
end

%------------- END OF CODE --------------