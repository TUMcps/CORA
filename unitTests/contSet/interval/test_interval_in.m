function res = test_interval_in
% test_interval_in - unit test function of in
%
% Syntax:  
%    res = test_interval_in
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

% Author:       Mark Wetzlinger, Adrian Kulmburg
% Written:      27-Sep-2019
% Last update:  12-March-2021 (MW, add empty case)
%               01-July-2021 (AK, merged test_interval_containsPoint with test_interval_in)
% Last revision:---

%------------- BEGIN CODE --------------

% 1. empty case
I = interval();
res_empty = in(I,I);

% 2. non-empty case
I = interval([-3;-2],[5;4]);
Z_in = zonotope([0.5, 2, 1;
                 0,   1,-0.7]);
Z_out = zonotope([6.5, 2, 1;
                 -3,   1,-0.7]);

res_in = in(I, Z_in);
res_out = in(I, Z_out); % should be out, but interval overapprox is in

res = res_empty && res_in && ~res_out;

% 3. point-containment case
% create interval
lower = [-3; -9; -4; -7; -1];
upper = [4;   2;  6;  3;  8];
Int = interval(lower, upper);

% points inside
p_inside = [0, 0, 0, 0, 0;
            1,-4, 3,-6, 5;
           -2,-6,-2, 2, 7;
           -3, 2, 6,-7, 8]';
res_inside = in(Int,p_inside);

% points outside
p_outside = [5, 3, 7, 4, 9;
             1,-4,-5,-6, 5;
            -2, 3,-2, 6, 7;
            -3, 2, 6,-7,10]';
res_outside = ~in(Int,p_outside);

% check if all points are in interval
res_point = res_inside && res_outside;

res = res & res_point;

if res
    disp('test_in successful');
else
    disp('test_in failed');
end

%------------- END OF CODE --------------